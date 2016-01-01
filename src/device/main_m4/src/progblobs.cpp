//
// begin license header
//
// This file is part of Pixy CMUcam5 or "Pixy" for short
//
// All Pixy source code is provided under the terms of the
// GNU General Public License v2 (http://www.gnu.org/licenses/gpl-2.0.html).
// Those wishing to use Pixy source code, software and/or
// technologies under different licensing terms should contact us at
// cmucam@cs.cmu.edu. Such licensing terms are available for
// all portions of the Pixy codebase presented here.
//
// end license header
//

#include "progblobs.h"
#include "pixy_init.h"
#include "camera.h"
#include "led.h"
#include "conncomp.h"
#include "serial.h"
#include "rcservo.h"
#include "exec.h"


bool g_ledSet = false;
static uint8_t g_state=0;

Program g_progBlobs =
{
	"Color_connected_components",
	"perform color connected components",
	blobsSetup, 
	blobsLoop
};


int blobsSetup()
{
	uint8_t c;

	// setup camera mode
	cam_setMode(CAM_MODE1);

	// if there have been any parameter changes, we should regenerate the LUT (do it regardless)
	g_blobs->m_clut.generateLUT();	
			
	// setup qqueue and M0
	g_qqueue->flush();
	exec_runM0(0);

	// flush serial receive queue
	while(ser_getSerial()->receive(&c, 1));

	g_state = 0; // reset recv state machine
	return 0;
}

void handleRecv()
{
	uint8_t i, a;
	static uint16_t w=0xffff;
	static uint8_t lastByte;
	uint16_t s0, s1;
	Iserial *serial = ser_getSerial();

	for (i=0; i<10; i++)
	{
		switch(g_state)
		{	
		case 0: // reset 
			lastByte = 0xff;  // This is not part of any of the sync word most significant bytes
			g_state = 1;
		 	break;

		case 1:	// sync word
			if(serial->receive(&a, 1))
			{
				w = lastByte << 8;
				w |= a;
				lastByte = a;
				g_state = 2;	// compare
			}
			break;

		case 2:	 // receive data byte(s)
			if (w==SYNC_SERVO)
			{	// read rest of data
				if (serial->receiveLen()>=4)
				{
					serial->receive((uint8_t *)&s0, 2);
					serial->receive((uint8_t *)&s1, 2);

					//cprintf("servo %d %d\n", s0, s1);
					rcs_setPos(0, s0);
					rcs_setPos(1, s1);

					g_state = 0;
				}
			}
			else if (w==SYNC_CAM_BRIGHTNESS)
			{
				if(serial->receive(&a, 1))
				{
					cam_setBrightness(a);
					g_state = 0;
				}
			}
			else if (w==SYNC_SET_LED)
			{
				if (serial->receiveLen()>=3)
				{
					uint8_t r, g, b;
					serial->receive(&r, 1);
					serial->receive(&g, 1);
					serial->receive(&b, 1);

					led_setRGB(r, g, b);
					//cprintf("%x %x %x\n", r, g ,b);

					g_ledSet = true; // it will stay true until the next power cycle
					g_state = 0;
				}
			}
			else 
				g_state = 1; // try another word, but read only a byte
			break;

		default:
			g_state = 0; // try another whole word
			break;
		}
	}
}

int blobsLoop()
{
#if 1
	BlobA *blobs;
	BlobB *ccBlobs;
	uint32_t numBlobs, numCCBlobs;
	static uint32_t drop = 0;

	// create blobs
	++drop;
	// give slow controller a chance to read a few more color objs by skipping blobification
	// the static count down var below is decremented further down in this function
	static int16_t skipFrameCntDwn = 0;
	if( skipFrameCntDwn<=0){
		if (g_blobs->blobify()<0)
		{
			DBG("drop %d\n", drop);
			return 0;
		}
	}else{
		// just read out the qqueue data from the M0
		// and update the data for the auto white and bright CCLs.
		// Don't setup any blobs.
		g_blobs->runlengthAnalysis(false);
	}
	// handle received data immediately
	if (g_interface!=SER_INTERFACE_LEGO)
		handleRecv();

	// reset AEC dead band to large default after adjusting AEC brightness set point
	if(g_brightCntDwn>0){
		if( (g_brightCntDwn-=1)<=0) cam_stabilizeBrightness();
	}

	static int16_t skipCnt=1;
	// run CCLs @ 25Hz as in pixymon
	if(skipCnt--<=0){
		skipCnt=1;

		// trigger auto brightness control loop
		uint8_t brght = 0;
		if(g_blobs->m_autoBrightGain>0.0f){
			brght = g_blobs->updateAutoBright();
			cam_setBrightness(brght); // does not access camera registers if brightness value has not changed since last call
		}

		// trigger auto white balance control loop
		static uint32_t wbvOld=0;
		uint32_t wbv=0;
		if(!cam_getAWB()){
			if(g_blobs->m_autoWhiteGain>0.0f) wbv=g_blobs->updateAutoWhite();
			if(wbv!=wbvOld){
				wbvOld=wbv;
				cam_setWBV(wbv);
			}
		}else{
			wbv = cam_getWBV();
			g_blobs->setAutoWhiteWBV(wbv);
			wbvOld=wbv;
		}

		/*static int16_t cnt=0;
		if(++cnt>10){
			cnt=0;
			DBG("brght=%d wbv=0x%06X\n", brght, wbv);
		}*/
	}

	// give slow controller a chance to read a few more color objs
	if( skipFrameCntDwn--<=0){
		skipFrameCntDwn = g_skipFrames;
		// send blobs
		g_blobs->getBlobs(&blobs, &numBlobs, &ccBlobs, &numCCBlobs);
		cc_sendBlobs(g_chirpUsb, blobs, numBlobs, ccBlobs, numCCBlobs);

		ser_getSerial()->update();
	}
	// if user isn't controlling LED, set it here, according to biggest detected object
	if (!g_ledSet)
		cc_setLED();
	
	// deal with any latent received data until the next frame comes in
	while(!g_qqueue->queued())
	{
		if (g_interface!=SER_INTERFACE_LEGO)
			handleRecv();
	}

#endif
#if 0
	Qval qval;
	int j = 0;
	static int i = 0;
	while(1)
	{
		if (g_qqueue->dequeue(&qval))
		{
			j++;
			if (qval.m_col>=0xfffe)
			{
				cprintf("%d: %d %x\n", i++, j, qval.m_col);
				break;
			}
		}
	}
#endif
#if 0
	BlobA *blobs;
	BlobB *ccBlobs;
	uint32_t numBlobs, numCCBlobs;
	static uint32_t drop = 0;

	// create blobs
	if (g_blobs->blobify()<0)
	{
		DBG("drop %d\n", drop++);
		return 0;
	}
	g_blobs->getBlobs(&blobs, &numBlobs, &ccBlobs, &numCCBlobs);
	cc_sendBlobs(g_chirpUsb, blobs, numBlobs, ccBlobs, numCCBlobs);

#endif

 	return 0;
}
