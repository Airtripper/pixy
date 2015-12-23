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

#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef PIXY
 #include "pixy_init.h"
 #include "misc.h"
 #include "debug.h"
#else
 #include "debug.h" // this and the one above are not the same!
#endif

#include "colorlut.h"
#include "calc.h"

#ifndef PIXY
bool g_logExp = true;
#endif

IterPixel::IterPixel(const Frame8 &frame, const RectA &region)
{
    m_frame = frame;
    m_region = region;
    m_points = NULL;
    reset();
}

IterPixel::IterPixel(const Frame8 &frame, const Points *points)
{
    m_frame = frame;
    m_points = points;
    reset();
}

bool IterPixel::reset(bool cleari)
{
    if (cleari)
        m_i = 0;
    if (m_points)
    {
        if (m_points->size()>m_i)
        {
            m_region = RectA((*m_points)[m_i].m_x, (*m_points)[m_i].m_y, CL_GROW_INC, CL_GROW_INC);
            m_i++;
        }
        else
            return false; // empty!
    }
    m_x = m_y = 0;
    m_pixels = m_frame.m_pixels + (m_region.m_yOffset | 1)*m_frame.m_width + (m_region.m_xOffset | 1);
    return true;
}

bool IterPixel::next(UVPixel *uv, RGBPixel *rgb, bool omitCutOnY)
{
    if (m_points)
    {
        if (nextHelper(uv, rgb, omitCutOnY))
            return true; // working on the current block
        else // get new block
        {
            if (reset(false)) // reset indexes, increment m_i, get new block
                return nextHelper(uv, rgb, omitCutOnY);  // we have another block!
            else
                return false; // blocks are empty
        }
    }
    else
        return nextHelper(uv, rgb, omitCutOnY);
}


bool IterPixel::nextHelper(UVPixel *uv, RGBPixel *rgb, bool omitCutOnY)
{
    int32_t r, g1, g2, b, u, v, c, miny=CL_MIN_Y;

    while(1)
    {
        if (m_x>=m_region.m_width)
        {
            m_x = 0;
            m_y += 2;
            m_pixels += m_frame.m_width*2;
        }
        if (m_y>=m_region.m_height)
            return false;

        r = m_pixels[m_x];
        g1 = m_pixels[m_x - 1];
        g2 = m_pixels[-m_frame.m_width + m_x];
        b = m_pixels[-m_frame.m_width + m_x - 1];
        if (rgb)
        {
            rgb->m_r = r;
            rgb->m_g = (g1+g2)/2;
            rgb->m_b = b;
        }
        if (uv)
        {
            c = r+g1+b;
            if (c<miny && !omitCutOnY)
            {
                m_x += 2;
                continue;
            }
            u = ((r-g1)<<CL_LUT_ENTRY_SCALE)/c;
            c = r+g2+b;
            if (c<miny && !omitCutOnY)
            {
                m_x += 2;
                continue;
            }
            v = ((b-g2)<<CL_LUT_ENTRY_SCALE)/c;

            uv->m_u = u;
            uv->m_v = v;
        }

        m_x += 2;
        return true;
    }
}

uint32_t IterPixel::averageRgb(uint32_t *pixels)
{
	RGBPixel rgb;
	uint32_t r, g, b, n;
	reset();
	for (r=g=b=n=0; next(NULL, &rgb); n++)
	{
		r += rgb.m_r;
		g += rgb.m_g;
		b += rgb.m_b;		
	}

	r /= n;
	g /= n;
	b /= n;

	if (pixels)
		*pixels = n;
	return (r<<16) | (g<<8) | b;
}

ColorLUT::ColorLUT(uint8_t *lut)
{
	int i; 
    m_lut = lut;
    m_useExpSigs = false;
    memset((void *)m_signatures, 0, sizeof(ColorSignature)*CL_NUM_SIGNATURES);
    memset((void *)m_runtimeSigs, 0, sizeof(RuntimeSignature)*CL_NUM_SIGNATURES);
    for(int i=0; i<CL_NUM_SIGNATURES; ++i) m_expSigs[i] = ExperimentalSignature();
    clearLUT();

    setMinBrightness(CL_DEFAULT_MINY);
    m_minRatio = CL_MIN_RATIO;
    m_maxDist = CL_MAX_DIST;
    m_ratio = CL_DEFAULT_TOL;
    m_ccGain = CL_DEFAULT_CCGAIN;
	for (i=0; i<CL_NUM_SIGNATURES; i++)
		m_sigRanges[i] = CL_DEFAULT_SIG_RANGE;

#ifndef PIXY
    m_useExpLut = false;
    ColorLutCalculatorExp::setupDivLut();
#endif
}


ColorLUT::~ColorLUT()
{
}

#if 0
void ColorLUT::calcRatios(IterPixel *ip, ColorSignature *sig, float ratios[])
{
    bool ubounded, vbounded;
    UVPixel uv;
    uint32_t un=0, vn=0, n=0, counts[4];
    longlong usum=0, vsum=0;
    counts[0] = counts[1] = counts[2] = counts[3] = 0;

    ip->reset();
    while(ip->next(&uv))
    {
        ubounded = true;
        vbounded = true;

        if (uv.m_u>sig->m_uMin)
            counts[0]++;
		else
            ubounded = false;

        if (uv.m_u<sig->m_uMax)
            counts[1]++;
		else
            ubounded = false;

        if (uv.m_v>sig->m_vMin)
            counts[2]++;
		else
            vbounded = false;

        if (uv.m_v<sig->m_vMax)
            counts[3]++;
		else
            vbounded = false;

        // only use pixels that are within our test bounds to form the mean
        if (ubounded)
        {
            usum += uv.m_u;
			un++;
		}
		if (vbounded)
		{
            vsum += uv.m_v;
            vn++;
        }
		n++;
    }

    // calc ratios
    ratios[0] = (float)counts[0]/n;
    ratios[1] = (float)counts[1]/n;
    ratios[2] = (float)counts[2]/n;
    ratios[3] = (float)counts[3]/n;
   // calc mean (because it's cheap to do it here)
    sig->m_uMean = usum/un;
    sig->m_vMean = vsum/vn;
 //	printf("%d %d %d %d %d %d %d %d %d\n", un, vn, n, sig->m_uMin, sig->m_uMax, sig->m_vMin, sig->m_vMax, sig->m_uMean, sig->m_vMean);
}

#else

void ColorLUT::calcRatios(IterPixel *ip, ColorSignature *sig, float ratios[])
{
    UVPixel uv;
    uint32_t n=0, counts[4];
    counts[0] = counts[1] = counts[2] = counts[3] = 0;

    ip->reset();
    while(ip->next(&uv))
    {
        if (uv.m_u>sig->m_uMin)
            counts[0]++;

        if (uv.m_u<sig->m_uMax)
            counts[1]++;

        if (uv.m_v>sig->m_vMin)
            counts[2]++;

        if (uv.m_v<sig->m_vMax)
            counts[3]++;

        n++;
    }

    // calc ratios
    ratios[0] = (float)counts[0]/n;
    ratios[1] = (float)counts[1]/n;
    ratios[2] = (float)counts[2]/n;
    ratios[3] = (float)counts[3]/n;
    // calc mean (because it's cheap to do it here)
    // the mean would be quite touchy to outliers
    // thus take the middle of the 80% limits.
    sig->m_uMean = (sig->m_uMin + sig->m_uMax)/2;
    sig->m_vMean = (sig->m_vMin + sig->m_vMax)/2;
}
#endif

void ColorLUT::iterate(IterPixel *ip, ColorSignature *sig)
{
    int32_t scale;
    float ratios[4];

    // binary search -- this rouine is guaranteed to find the right value +/- 1, which is good enough!
    // find all four values, umin, umax, vmin, vmax simultaneously
    for (scale=1<<30, sig->m_uMin=sig->m_uMax=sig->m_vMin=sig->m_vMax=0; scale!=0; scale>>=1)
    {
        calcRatios(ip, sig, ratios);
        if (ratios[0]>m_ratio)
            sig->m_uMin += scale;
        else
            sig->m_uMin -= scale;

        if (ratios[1]>m_ratio)
            sig->m_uMax -= scale;
        else
            sig->m_uMax += scale;

        if (ratios[2]>m_ratio)
            sig->m_vMin += scale;
        else
            sig->m_vMin -= scale;

        if (ratios[3]>m_ratio)
            sig->m_vMax -= scale;
        else
            sig->m_vMax += scale;
    }
}




int ColorLUT::generateSignature(const Frame8 &frame, const RectA &region, uint8_t signum)
{
    if (signum<1 || signum>CL_NUM_SIGNATURES)
        return -1;
    // this is cool-- this routine doesn't allocate any extra memory other than some stack variables
    IterPixel ip(frame, region);
    iterate(&ip, m_signatures+signum-1);
    m_signatures[signum-1].m_type = 0;

    ip.reset();
    accExpSig(signum).init( ip);

    updateSignature(signum);
    return 0;
}


int ColorLUT::generateSignature(const Frame8 &frame, const Point16 &point, Points *points, uint8_t signum)
{
	if (signum<1 || signum>CL_NUM_SIGNATURES)
		return -1;
    // this routine requires some memory to store the region which consists of some consistently-sized blocks
    growRegion(frame, point, points);
    IterPixel ip(frame, points);
    iterate(&ip, m_signatures+signum-1);
	m_signatures[signum-1].m_type = 0;

    ip.reset();
    accExpSig(signum).init( ip);

    updateSignature(signum);
    return 0;
}

void ColorLUT::updateSignature(uint8_t signum)
{
    float range;

	if (signum<1 || signum>CL_NUM_SIGNATURES)
		return;
	signum--;

    if (m_signatures[signum].m_type==CL_MODEL_TYPE_COLORCODE)
        range = m_sigRanges[signum]*m_ccGain;
	else
		range = m_sigRanges[signum];
    m_runtimeSigs[signum].m_uMin = m_signatures[signum].m_uMean + (m_signatures[signum].m_uMin - m_signatures[signum].m_uMean)*range;
	m_runtimeSigs[signum].m_uMax = m_signatures[signum].m_uMean + (m_signatures[signum].m_uMax - m_signatures[signum].m_uMean)*range;
	m_runtimeSigs[signum].m_vMin = m_signatures[signum].m_vMean + (m_signatures[signum].m_vMin - m_signatures[signum].m_vMean)*range;
	m_runtimeSigs[signum].m_vMax = m_signatures[signum].m_vMean + (m_signatures[signum].m_vMax - m_signatures[signum].m_vMean)*range;

    m_runtimeSigs[signum].m_rgbSat = saturate(m_signatures[signum].m_rgb);
}

ColorSignature *ColorLUT::getSignature(uint8_t signum)
{
	if (signum<1 || signum>CL_NUM_SIGNATURES)
		return NULL;

	return m_signatures+signum-1;
}

int ColorLUT::setSignature(uint8_t signum, const ColorSignature &sig)
{
	if (signum<1 || signum>CL_NUM_SIGNATURES)
		return -1;

	m_signatures[signum-1] = sig;
	updateSignature(signum);
	return 0;
}


int ColorLUT::generateLUT()
{

#ifdef PIXY
    uint32_t timer;
    setTimer(&timer);
#if 1
    const uint32_t keepAliveTmOut = 50000;
    uint32_t keepAliveTmr;
    setTimer(&timer);
#endif
#else
    int collisions = 0;
#endif

    clearLUT();

    if(m_useExpSigs){

        float gValMinF=23.f;
        for (uint16_t s=1; s<=CL_NUM_SIGNATURES; ++s){
            const ExperimentalSignature& es = expSig(s);
            if(es.isActive() && es.hsvValMin()<gValMinF) gValMinF=es.hsvValMin();
        }
        int16_t gValMin = (gValMinF*255.f)+0.5f;

        // loop with sufficient granularity thru uv space
        const int16_t bm=(1<<8)-1;
        const int16_t stpG=3;
#ifdef PIXY
        const int16_t off=3;
        const int16_t stpUV=8;
#else
        int16_t off = m_useExpLut ? 0 : 3;
        int16_t stpUV=m_useExpLut ? 3 : 8;
#endif

        for( int16_t u=-bm+off; u<=bm; u+=stpUV){
            for( int16_t v=-bm+off; v<=bm; v+=stpUV){

                // prepare limits for the inner loop thru g values
                // 0 <= {r=u+g,b=v+g,g} <=255
                int16_t gMin=0;
                if(gMin<-u)gMin=-u;
                if(gMin<-v)gMin=-v;
                int16_t gMax=bm;
                if(gMax>bm-u)gMax=bm-u;
                if(gMax>bm-v)gMax=bm-v;

                 // calc the (u,v) position in the color LUT
                int16_t ui = u >> (9-CL_LUT_COMPONENT_SCALE);
                ui &= (1<<CL_LUT_COMPONENT_SCALE)-1;
                int16_t vi = v >> (9-CL_LUT_COMPONENT_SCALE);
                vi &= (1<<CL_LUT_COMPONENT_SCALE)-1;
                int16_t lutIdx = (ui<<CL_LUT_COMPONENT_SCALE)+ vi;

                for( uint16_t g=gMin; g<=gMax; g+=stpG){
#ifdef PIXY
                    if(getTimer(keepAliveTmr)>keepAliveTmOut){
                         g_chirpUsb->service(); // keep alive
                         setTimer(&keepAliveTmr);
                    }
#endif
                    // bail out early on dark colours (guess the speed gain is small, but now it's in)
                    if(g<gValMin && u+g<gValMin && v+g<gValMin)continue;

                    // determine the most probable signature by comparing distances in the u/v plane
                    float gf=g*rgbNorm;
                    float rf=(u+g)*rgbNorm;
                    float bf=(v+g)*rgbNorm;
                    float dst2Min = 23.0;
                    uint8_t bestSigId = 0;
                    for (uint16_t s=1; s<=CL_NUM_SIGNATURES; ++s){
                        // check signature compatibility and calc the distance in the (u,v) plane
                        float uf,vf;
                        const ExperimentalSignature& sig = expSig(s);
                        if( sig.isActive() && sig.isRgbAccepted(rf,gf,bf, uf,vf)){
                            float du = uf-sig.uMed();
                            float dv = vf-sig.vMed();
                            float dst2 = du*du+dv*dv;
                            if(dst2<dst2Min){
                                dst2Min=dst2;
                                bestSigId=s;
                            }
                        }
                    }
                    if(bestSigId){
#ifndef PIXY
                        if(m_lut[lutIdx] & ~(1<<(bestSigId-1))) ++collisions;

                        if( m_useExpLut ){
                            int16_t ui,vi;
                            ColorLutCalculatorExp::calcUV( u+g,g,v+g, ui,vi);
                            // ... and store a 7bit bitmap of compatible signatures in there
                            // Usually not more than one bit should should be set,
                            // but nearby signatures might overlap due to the non perfect
                            // division approximation used in the M0 preselection (u,v) caculation.
                            // Those collisions are re-checked in the final filter step performed on the M4
                            uint16_t lutIdx = (vi<<CL_LUT_COMPONENT_SCALE) | ui; // alternative LUT arrangement that can be visualized, see ColorLutCalculatorExp::calcUV
                            m_lut[lutIdx] |= 1<<(bestSigId-1);
                        }else
#endif
                        {
                            // update the 7bit bitmap of compatible signatures in the uv-LUT
                            m_lut[lutIdx] |= 1<<(bestSigId-1);

                            // Update the y-LUT
                            // store for each value of the 7-signature-bitmap ([1..127]) entered in the uv-LUT
                            // the minimum and maximum y=r+g+b=u+v+3g seen during the LUT setup sweep
                            // As Pixy accumulates two pixels per qVal, multiply the min and max with 2.
                            uint32_t y = 2*(u+v+3*g);
                            uint16_t* yLut = (uint16_t*)(m_lut+0x1000);
                            uint32_t yLutIdx = 2*m_lut[lutIdx];
                            if(y<yLut[yLutIdx]) yLut[yLutIdx]=y;
                            if(y>yLut[yLutIdx+1]) yLut[yLutIdx+1]=y;
                        }
                    }
                }
            }
        }
    }
    else
    {
        int32_t r, g, b, u, v, y, bin, sig;

        // recalc bounds for each signature
        for (r=0; r<CL_NUM_SIGNATURES; r++)
            updateSignature(r+1);

        for (r=0; r<1<<8; r+=1<<(8-CL_LUT_COMPONENT_SCALE))
        {
            for (g=0; g<1<<8; g+=1<<(8-CL_LUT_COMPONENT_SCALE))
            {
                for (b=0; b<1<<8; b+=1<<(8-CL_LUT_COMPONENT_SCALE))
                {
                    y = r+g+b;

                    if (y<(int32_t)m_miny)
                        continue;
                    u = ((r-g)<<CL_LUT_ENTRY_SCALE)/y;
                    v = ((b-g)<<CL_LUT_ENTRY_SCALE)/y;

                    for (sig=0; sig<CL_NUM_SIGNATURES; sig++)
                    {
                        if (m_signatures[sig].m_uMin==0 && m_signatures[sig].m_uMax==0)
                            continue;
                        if ((m_runtimeSigs[sig].m_uMin<u) && (u<m_runtimeSigs[sig].m_uMax) &&
                                (m_runtimeSigs[sig].m_vMin<v) && (v<m_runtimeSigs[sig].m_vMax))
                        {
                            int32_t u = r-g; // evillive see below (made local to get a reliable collision count)
                            u >>= 9-CL_LUT_COMPONENT_SCALE;
                            u &= (1<<CL_LUT_COMPONENT_SCALE)-1;
                            int32_t v = b-g; // evillive see below (made local to get a reliable collision count)
                            v >>= 9-CL_LUT_COMPONENT_SCALE;
                            v &= (1<<CL_LUT_COMPONENT_SCALE)-1;

                            bin = (u<<CL_LUT_COMPONENT_SCALE)+ v;
#ifndef PIXY
                            if(m_lut[bin] && m_lut[bin]!= (1<<sig) ) ++collisions;
#endif

                            if (m_lut[bin]==0 || m_lut[bin] > (1<<sig) ){
                                // lower index signatures have higher prio and kick lower prio signatures out of the LUT
                                m_lut[bin] = 1<<sig;
#ifdef PIXY // don't bail out in cooked mode as it would affect the collision counting
                                break; // bail out here as higher index signatures don't have a chance to get into this LUT entry any more
#endif
                            }
                            // Overwriting u and v and not bailing out of the loop here
                            // is either a bug or an ELO 2000+ design i don't understand.
                            // Not absolutely sure, but I think it doesn't harm.
                        }
                    }
                }
            }
        }
    }

#ifndef PIXY
    EXPLOG("LUT Dump (collisions=%d)", collisions);

    uint16_t* yLut = (uint16_t*)(m_lut+0x1000);
    for(uint32_t i=0; i<128; ++i) if(yLut[2*i]!=0xffff && yLut[2*i+1]!=0xffff) EXPLOG("7sigs=%2X  %3d < y < %3d", i, yLut[2*i], yLut[2*i+1] );

    const int sz = (1<<CL_LUT_COMPONENT_SCALE);
    const int strLen = sz*4+1;
    char str[strLen];
    for(int v=sz-1; v>=0; --v){
        unsigned int pos=0;
        for(int u=0; u<sz; ++u){
            int ut,vt;
            if(m_useExpLut){
                ut=u;
                vt=v;
            }else{
                ut = (v-(1<<(CL_LUT_COMPONENT_SCALE-1))) & ((1<<CL_LUT_COMPONENT_SCALE)-1);
                vt = (u-(1<<(CL_LUT_COMPONENT_SCALE-1))) & ((1<<CL_LUT_COMPONENT_SCALE)-1);
            }
            if(m_lut[vt*sz+ut])
                pos += snprintf( str+pos, (pos<strLen ? strLen-pos : 0), "%2X", m_lut[vt*sz+ut]);
            else
                pos += snprintf( str+pos, (pos<strLen ? strLen-pos : 0), "..");
        }
        EXPLOG(str);
    }
#endif

#ifdef PIXY
    DBG("genLUT expSig=%d %dms", m_useExpSigs, (getTimer(timer)+500)/1000);
#endif

    return 0;
}


void ColorLUT::clearLUT(uint8_t signum)
{
    int i;

    for (i=0; i<CL_LUT_SIZE; i++)
    {
        if (signum==0)
            m_lut[i] = 0;
        else
            m_lut[i] &= ~(1<<(signum-1));
    }

    // dirty !
    // reset the yLUT
    uint16_t* yLut = (uint16_t*)(m_lut+0x1000);
    for(uint32_t i=0; i<128; ++i){
        yLut[2*i]= m_useExpSigs ? 0xffff : 0;
        yLut[2*i+1]= m_useExpSigs ? 0 : 0xffff;
    }
}


bool ColorLUT::growRegion(RectA *region, const Frame8 &frame, uint8_t dir)
{
    if (dir==0) // grow left
    {
        if (region->m_xOffset>=CL_GROW_INC)
        {
            region->m_xOffset -= CL_GROW_INC;
            region->m_width += CL_GROW_INC;
        }
        else
            return true;
    }
    else if (dir==1) // grow top
    {
        if (region->m_yOffset>=CL_GROW_INC)
        {
            region->m_yOffset -= CL_GROW_INC;
            region->m_height += CL_GROW_INC;
        }
        else
            return true;
    }
    else if (dir==2) // grow right
    {
        if (region->m_xOffset+region->m_width+CL_GROW_INC>frame.m_width)
            return true;
        region->m_width += CL_GROW_INC;
    }
    else if (dir==3) // grow bottom
    {
        if (region->m_yOffset+region->m_height+CL_GROW_INC>frame.m_height)
            return true;
        region->m_height += CL_GROW_INC;
    }
    return false;
}


float ColorLUT::testRegion(const RectA &region, const Frame8 &frame, UVPixel *mean, Points *points)
{
    UVPixel subMean;
    float distance;
    RectA subRegion(0, 0, CL_GROW_INC, CL_GROW_INC);
    subRegion.m_xOffset = region.m_xOffset;
    subRegion.m_yOffset = region.m_yOffset;
    bool horiz = region.m_width>region.m_height;
    uint32_t i, test, endpoint = horiz ? region.m_width : region.m_height;

    for (i=0, test=0; i<endpoint; i+=CL_GROW_INC)
    {
        getMean(subRegion, frame, &subMean);
#ifdef PIXY
    distance = vsqrtf((float)((mean->m_u-subMean.m_u)*(mean->m_u-subMean.m_u) + (mean->m_v-subMean.m_v)*(mean->m_v-subMean.m_v)));
#else
    distance = sqrtf((float)((mean->m_u-subMean.m_u)*(mean->m_u-subMean.m_u) + (mean->m_v-subMean.m_v)*(mean->m_v-subMean.m_v)));
#endif
        if ((uint32_t)distance<m_maxDist)
        {
            int32_t n = points->size();
            mean->m_u = ((longlong)mean->m_u*n + subMean.m_u)/(n+1);
            mean->m_v = ((longlong)mean->m_v*n + subMean.m_v)/(n+1);
            if (points->push_back(Point16(subRegion.m_xOffset, subRegion.m_yOffset))<0)
                break;
            //DBG("add %d %d %d", subRegion.m_xOffset, subRegion.m_yOffset, points->size());
            test++;
        }

        if (horiz)
            subRegion.m_xOffset += CL_GROW_INC;
        else
            subRegion.m_yOffset += CL_GROW_INC;
    }

    //DBG("return %f", (float)test*CL_GROW_INC/endpoint);
    return (float)test*CL_GROW_INC/endpoint;
}


void ColorLUT::growRegion(const Frame8 &frame, const Point16 &seed, Points *points)
{
    uint8_t dir, done;
    RectA region, newRegion;
    UVPixel mean;
    float ratio;

    done = 0;

    // create seed 2*CL_GROW_INCx2*CL_GROW_INC region from seed position, make sure it's within the frame
    region.m_xOffset = seed.m_x;
    region.m_yOffset = seed.m_y;
    if (growRegion(&region, frame, 0))
        done |= 1<<0;
    else
        points->push_back(Point16(region.m_xOffset, region.m_yOffset));
    if (growRegion(&region, frame, 1))
        done |= 1<<1;
    else
        points->push_back(Point16(region.m_xOffset, region.m_yOffset));
    if (growRegion(&region, frame, 2))
        done |= 1<<2;
    else
        points->push_back(Point16(seed.m_x, region.m_yOffset));
    if (growRegion(&region, frame, 3))
        done |= 1<<3;
    else
        points->push_back(seed);

    getMean(region, frame, &mean);

    while(done!=0x0f)
    {
        for (dir=0; dir<4; dir++)
        {
            newRegion = region;
            if (done&(1<<dir))
                continue;
            else if (dir==0) // left
                newRegion.m_width = 0;
            else if (dir==1) // top
                newRegion.m_height = 0; // top and bottom
            else if (dir==2) // right
            {
                newRegion.m_xOffset += newRegion.m_width;
                newRegion.m_width = 0;
            }
            else if (dir==3) // bottom
            {
                newRegion.m_yOffset += newRegion.m_height;
                newRegion.m_height = 0;
            }

            if (growRegion(&newRegion, frame, dir))
                done |= 1<<dir;
            else
            {
                ratio = testRegion(newRegion, frame, &mean, points);
                if (ratio<m_minRatio)
                    done |= 1<<dir;
                else
                    growRegion(&region, frame, dir);
            }
        }
    }
}


void ColorLUT::getMean(const RectA &region ,const Frame8 &frame, UVPixel *mean)
{
    UVPixel uv;
    uint32_t n=0;
    IterPixel ip(frame, region);

    longlong usum=0, vsum=0;

    while(ip.next(&uv))
    {
        usum += uv.m_u;
        vsum += uv.m_v;
        n++;
    }

    mean->m_u = usum/n;
    mean->m_v = vsum/n;
}

void ColorLUT::setSigRange(uint8_t signum, float range)
{
	if (signum<1 || signum>CL_NUM_SIGNATURES)
		return;
	m_sigRanges[signum-1] = range;
}

void ColorLUT::setGrowDist(uint32_t dist)
{
	m_maxDist = dist;
}

void ColorLUT::setMinBrightness(float miny)
{
    m_miny = 3*((1<<8)-1)*miny;
    if (m_miny==0)
        m_miny = 1;

}

void ColorLUT::setCCGain(float gain)
{
    m_ccGain = gain;
}

uint32_t ColorLUT::getType(uint8_t signum)
{
    if (signum<1 || signum>CL_NUM_SIGNATURES)
        return 0;

    return m_signatures[signum-1].m_type;
}

#if 0
uint32_t ColorLUT::getColor(uint8_t signum)
{
    int32_t r, g, b, max, u, v;

    if (signum<1 || signum>CL_NUM_SIGNATURES)
        return 0;

    u = m_signatures[signum-1].m_uMean;
    v = m_signatures[signum-1].m_vMean;

    // u = r-g
    // v = b-g
    if (abs(u)>abs(v))
    {
        if (u>0)
        {
            r = u;
            if (v>0)
                g = 0;
            else
                g = -v;
            b = v+g;
        }
        else
        {
            g = -u;
            r = 0;
            b = v+g;
        }
    }
    else
    {
        if (v>0)
        {
            b = v;
            if (u>0)
                g = 0;
            else
                g = -u;
            r = u+g;
        }
        else
        {
            g = -v;
            b = 0;
            r = u+g;
        }
    }

    if (r>g)
        max = r;
    else
        max = g;
    if (b>max)
        max = b;

    // normalize
    if (max>0)
    {
        r = (float)r/max*255;
        g = (float)g/max*255;
        b = (float)b/max*255;
        return (r<<16) | (g<<8) | b;
    }
    else
        return 0;
}
#endif



#ifndef PIXY
uint8_t ColorLutCalculatorExp::s_divLut[256];

void ColorLutCalculatorExp::calcUV(int16_t r, int16_t g, int16_t b, int16_t &u, int16_t &v)
{
    // Approximate classic YUV luminance Y and boost to maximum range that fits into int16
    // 2:4:1 isn't a that bad approximation of 0.299:0.587:0.114
    uint8_t rs=5;
    uint8_t gs=6;
    uint8_t bs=4;
    int16_t y = (r<<rs) + (g<<gs) + (b<<bs);

    // Extend the r and b ranges accordingly
    // Do similar shifting as above as the LPC43xx's M0 doesn't support fast multiplication
    int16_t re = (r<<rs) + (r<<gs) + (r<<bs);  // r16*(32+64+16)=r16*112
    int16_t be = (b<<rs) + (b<<gs) + (b<<bs);

    // Calculate classic YUV chrominance values	u and v
    u = be-y;
    v = re-y;

    // tune the chroma vals to correct a bit the weight approximation errors above
    // and - more important - to cover as many bins of the LUT as possible
    u += (u>>3)+(u>>6);
    v += (v>>2)+(v>>4)+(v>>5)+(v>>6);

    // "Normalize" the chroma vals to the maximum of the rgb vals
    // This moves u/v vals of the gloomy but well saturated pixels from
    // the center to the bright edges of uv-plane, leaving only the low
    // saturated (gloomy and bright) in the center region.
    // first get the maximum of rgb ...
    short rgbm = r;
    if(g>rgbm) rgbm=g;
    if(b>rgbm) rgbm=b;

    // ... and approximate the division u/rgbm and v/rgbm
    // As division is even slower than multiplication
    // Thus do the following shift/add/shift approximation instead of a real division
    // shift values are taken from a look up table
    uint8_t lut = s_divLut[rgbm];
    uint8_t s1 = lut&0xf;
    if(s1){
        u += u>>s1;
        v += v>>s1;
    }
    uint8_t s2 = lut>>4;
    u >>= s2;
    v >>= s2;

    // first make them positive by a simple add
    // ... same number of cycles but less headache than the original shift method
    u += 1<<(CL_LUT_COMPONENT_SCALE-1);
    v += 1<<(CL_LUT_COMPONENT_SCALE-1);

    // a final range check that should always pass ... evillive could be removed (?)
    const int16_t rangeMask = (1<<CL_LUT_COMPONENT_SCALE)-1;
    if( u & ~rangeMask || v & ~rangeMask ) EXPLOG("Exp LUT fail %d %d %d %i %i", r,g,b,u,v);
    // but be defensive and limit u,v to the allowed range
    u &= rangeMask;
    v &= rangeMask;
}

void ColorLutCalculatorExp::setupDivLut()
{
    // Has been tweaked! Now it divides by b<<2
    // to get directly within the range of the color LUT
    const int16_t tweak = 8-CL_LUT_COMPONENT_SCALE;

    // A magic reference nominator seed for error estimation
    const int16_t nomi = 11614;
    s_divLut[0]=0xf0;  // x/0=0
    s_divLut[1]=(tweak<<4);

    float maxErr=0.0f;
    for(int16_t dnomi=2; dnomi<256; ++dnomi){ // loop through all denominators

        float refQ = float(nomi) / (dnomi<<tweak); // the "exact" result
        uint8_t msb = 1;
        while(dnomi>>msb) ++msb;	// pos of most significant denominator bit, used as loop limit below
        uint8_t s1MinErr = 0;       // stores first and second shift
        uint8_t s2MinErr = 0;       // for the
        float minErr = 23.0f;		// lowest quotient error found so far
        //int16_t qMinErr = 0; // just debug stuff

        // loop through reasonable bit shift distances
        // and search for the combination with the lowest quotient approximation error
        for(uint8_t s2=msb-1; s2<=msb+1+tweak; ++s2){
            for(uint8_t s1=0; s1<8; ++s1){

                int16_t apprxQ = nomi;
                if(s1) apprxQ += (nomi>>s1);
                apprxQ >>= s2;

                float err = fabs(apprxQ/refQ-1.0f);
                if(err<minErr){
                    minErr = err;
                    s1MinErr = s1;
                    s2MinErr = s2;
                    //qMinErr = apprxQ;
                }
            }
        }
        // store the best s1/s2 combination to the LUT
        s_divLut[dnomi]=s1MinErr |(s2MinErr<<4);

        //printf("%3u %2u %u %2x %5i %5.1f %2.1f%% %3.0f\n", dnomi, s1MinErr, s2MinErr, divLUT[dnomi], qMinErr, refQ, minErr*100.0f, qMinErr-refQ);
        if(minErr>maxErr)maxErr=minErr;
    }
    EXPLOG("Division LUT: seed=%u err=%.1f%%",nomi,maxErr*100.0f);
}
#endif
