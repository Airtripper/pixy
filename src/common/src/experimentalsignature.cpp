#include "experimentalsignature.h"
#include "colorlut.h"  // cyclic due to IterPixel

#include <math.h>

// parameter names or their format strings
const char* parName_eSigUse = "Use Experimental Signatures";
const char* parName_eSigPos = "eSigPos%d";
const char* parName_eSigAct = "eSig%d Active";
const char* parName_eSigHueRng = "eSig%d HueRange";
const char* parName_eSigSatMin = "eSig%d SatMin";
const char* parName_eSigSatMax = "eSig%d SatMax";
const char* parName_eSigValMin = "eSig%d ValMin";
const char* parName_eSigValMax = "eSig%d ValMax";


ExperimentalSignature::ExperimentalSignature():   
    m_posUV(),
    m_hsvSatMed(0.00f),  // the default is white
    m_hsvValMin(0.95f),
    m_hsvValMax(1.00f),
    m_hsvSatMin(0.00f),
    m_hsvSatMax(0.05f),
    m_hsvCosDeltaHueMin(1.0f),
    m_hsvHueRange(0.0f),
    m_isActive(false)
{}

ExperimentalSignature::~ExperimentalSignature()
{}


bool ExperimentalSignature::isRgbAccepted(float r, float g, float b, float& u, float& v) const
{
    // The colorspace used in here somehow trial and error evolved from the one
    // originally used in e.g. Blobs::runlengthAnalysis, IterPixel::nextHelper and ColorLUT::generateLUT of FW 2.0.8
    // I fact it finally applies a cut in the HSV color space
    // One might argue that calculating the HSV vector directly would be cheaper
    // good  point, maybe, but one wouldn't have the nice symmetric circle like HSV distribution
    // A rectangular distribution would make selections at the edges and corners difficult.
    // One should investigate which approach is better, but first check how the algorithm performs.
    // evillive: check efficiency and optimize!  up to three divisions and a squareroot  => 4*14 cycles + ? ~ 100 cycles?
    // utilize some of the M4 CMSIS-DSP stuff?
    // check: It looks like that the YUV (u,v)-distribution is nothing else than a projection of the rgb cube alongside the saturation=0 body diagonal.
    // Haven't found a confirmation of this yet. ... Use M4 CMSIS-DSP matrix/vector features for doing this?
    // this is more or less a duplicate of RuntimeSignature::calculateUV
    // clean this up, ... maybe at least some inline functions

    // calc rgb min and max
    float hsvVal = g; // aka rgbMax
    float rgbMin = g;
    if(r>hsvVal)hsvVal=r; else rgbMin=r;
    if(b>hsvVal)hsvVal=b;
    if(b<rgbMin)rgbMin=b;

    // might happen, floating point or g1/g2 bayer mess caused
//    if(hsvVal>1.0f){
//        hsvVal=1.0f;
//        //EXPLOG("val ouch");
//    }

    // bail out early if not within hsv value limits
    // and handling of black pixels (division by zero prevention)
    if( hsvVal<m_hsvValMin+bite || (hsvVal>m_hsvValMax && hsvVal<=1.0f) ){
        u=v=0.0f;
        return hsvVal+m_hsvValMin+m_hsvSatMin < bite;  // return true if signature accepts black
    }

    //float hsvValInv=1.0f/hsvVal;                     // save a division
    //float hsvSat = float(hsvVal-rgbMin) * hsvValInv; // save a division
    float hsvSat = float(hsvVal-rgbMin);               // save a division

    // might happen (?), floating point or g1/g2 bayer mess caused ... handled in the condition below
//    if(hsvSat>1.0f){
//        hsvSat=1.0f;
//        //EXPLOG("sat ouch");
//    }

    // bail out early if not within hsv saturation limits
    // and handling of grey pixels (division by zero prevention)
    //if( hsvSat<m_hsvSatMin+bite || (hsvSat>m_hsvSatMax && hsvSat<=1.0f) ) {
    if( hsvSat<m_hsvSatMin*hsvVal+bite || hsvSat>m_hsvSatMax*hsvVal) {
        u=v=0;
        return hsvSat+m_hsvSatMin<bite ; // return true if signature accepts 255 shades of grey ... :D
    }

    float hsvValInv=1.0f/hsvVal;    // save a division

    // calc classic YUV luminance y
    float y = yuv_wr*r + yuv_wg*g + yuv_wb*b;

    // calc classic YUV chrominance u and v and normalize them
    u = b-y;
    u *= hsvValInv * yuv_un;
    v = r-y;
    v *= hsvValInv * yuv_vn;

    // now we have a (u,v) vector pointing into a nice YUV hexagon in the u-v plane
    // HSV saturation grows linarly from center to the edges
    // u-v coordinates are independent of the HSV value
    // now make the hexagon a circle
    // we would need the sqrt anyway, but it adds another division
    float circleFac = hsvSat / sqrtf(u*u+v*v);  // evillive check if M4 FPU is used
    u *= circleFac;
    v *= circleFac;

    // We got here with with zero signature saturation => Accept the candidate
    // This bail out prevents a division by zero below.
    // We could have bailed out earlier, but we wanna return valid u and v values
    if(m_hsvSatMed<bite)return true;

    // calculate cosine delta hue using dot product: (u,v)_sig . (u,v)_pix / |(u,v)_sig| / |(u,v)_pix|
    float cosHue = u*uMed() + v*vMed();
    //cosHue /= hsvSat * m_hsvSatMed; // save one division

    //if(cosHue>0.99) EXPLOG("u=%.2f v=%.2f um=%.2f vm=%.2f c=%.24f ", u,v,m_uMed ,m_vMed,cosHue);

    // Are we within HSV hue limits?
    // return cosHue>m_hsvCosDeltaHueMin; // save one division
    return cosHue > m_hsvCosDeltaHueMin * hsvSat * m_hsvSatMed; // safe one division
}


void ExperimentalSignature::init( IterPixel& pixIter)
{
    // approximation of the hue and saturation median
    // do it with simple histograms
    // based on "Binapprox" algorithm described in
    // "Fast Computation of the Median by Successive Binning"
    // Ryan J. Tibshirani, Dept. of Statistics, Stanford University""

    // As calculating mean and sigma of the circular hue
    // might lead to surprising results, if it's distribution
    // peaks at the jump from +pi to -pi, we have to rotate
    // the hue peak into a save region
    // 1) abuse the histos, used later for hue and saturation, to calculate a mean u and v
    Histo hueHist( -1.0f, 1.0f);
    Histo satHist( -1.0f, 1.0f);
    RGBPixel rgbPix;
    pixIter.reset();
    while(pixIter.next( 0 ,&rgbPix)){
         float u, v, sat, val;
         translateRGB( rgbPix.m_r,rgbPix.m_g,rgbPix.m_b,  u, v, sat, val );
         hueHist.add( u);
         satHist.add( v);
    }
    float uMean = hueHist.mean();
    float vMean = satHist.mean();

    // calculate the polar angle of the mean (u,v) vector
    float hueOff = fabs(uMean)>bite && fabs(vMean)>bite ? atan2(vMean,uMean) : 0.0f;

    // now analyse the saturation and hue distribution
    // hue is rotated by hueOff to the center of the histo before
    hueHist.reset( -pi, pi);
    satHist.reset( 0.0f, 1.0f);
    pixIter.reset();
    while(pixIter.next( 0 ,&rgbPix)){
         float u, v, sat, val;
         translateRGB( rgbPix.m_r,rgbPix.m_g,rgbPix.m_b,  u, v, sat, val );
         float hue0 = sat>bite ?  atan2f(v,u) : 0.0f;
         hue0-=hueOff;
         if(hue0<-pi)hue0+=2.0f*pi;
         else if(hue0>pi)hue0-=2.0f*pi;
         hueHist.add( hue0);
         satHist.add( sat);
    }

    // the selection cuts on the value are initialized to fix defaults
    m_hsvValMin = 0.2f;
    m_hsvValMax = 1.0f;

    // saturation and hue selection cuts are defined by the values
    // for 5% and 95% of the cumulative distribution
    const float outLim = 0.05f;
    m_hsvSatMin = satHist.X( outLim)-outLim;
    if(m_hsvSatMin<0.0f) m_hsvSatMin=0.0f;
    //m_hsvSatMax = satHist.X( 1.0f-outLim);
    m_hsvSatMax = 1;
    m_hsvHueRange = 0.5f * ( hueHist.X(1.0f-outLim) - hueHist.X(outLim) );
    m_hsvCosDeltaHueMin = cosf( m_hsvHueRange );

    // now have a closer look at the hue and sat distributions
    // their median should be within a 1 one sigma window arround the mean val
    float mean = hueHist.mean();
    float sigma = hueHist.sigma();
    hueHist.reset( mean-sigma, mean+sigma);

    mean = satHist.mean();
    sigma = satHist.sigma();
    satHist.reset( mean-sigma, mean+sigma);

    pixIter.reset();
    while(pixIter.next( 0 ,&rgbPix)){
         float u, v, sat, val;
         translateRGB( rgbPix.m_r,rgbPix.m_g,rgbPix.m_b,  u, v, sat, val );
         float hue0 = sat>bite ?  atan2f(v,u) : 0.0f;
         hue0-=hueOff;
         if(hue0<-pi)hue0+=2.0f*pi;
         else if(hue0>pi)hue0-=2.0f*pi;
         hueHist.add( hue0);
         satHist.add( sat);
    }

    // now get the saturation median
    m_hsvSatMed=satHist.X(0.5);

    // and calculate the u, v and hue medians
    float hueMed=hueHist.X(0.5)+hueOff;
    m_posUV.m_uMed = m_hsvSatMed*cosf( hueMed);
    m_posUV.m_vMed = m_hsvSatMed*sinf( hueMed);
#ifndef PIXY
    EXPLOG("init: uMed=%.2f vMed=%.2f satMed=%.2f satMin=%.2f satMax=%.2f hueMed=%.2f hueRng=%.2f hueCos=%.4f hueOff=%.2f",
           uMed(), vMed(), m_hsvSatMed, m_hsvSatMin, m_hsvSatMax, hueMed*r2d, m_hsvHueRange*r2d, m_hsvCosDeltaHueMin, hueOff*r2d );
#endif
    m_isActive = true;
}

void ExperimentalSignature::translateRGB(float r, float g, float b, float &u, float &v, float& hsvSat, float& hsvVal)
{
    // this is more or less a duplicate of RuntimeSignature::isRgbAccepted
    // docu can be found there
    // clean this up!

    // calc rgb min and max
    hsvVal = g; // aka rbgMax
    float rgbMin = g;
    if(r>hsvVal)hsvVal=r; else rgbMin=r;
    if(b>hsvVal)hsvVal=b;
    if(b<rgbMin)rgbMin=b;

    // a black pixel => set u and v zero and bail out (division by zero prevention)
    if( hsvVal<bite){
        u = v = 0.0f;
        hsvSat = 0.0f;
        return;
    }

    float hsvValInv=1.0f/hsvVal;
    hsvSat = (hsvVal-rgbMin) * hsvValInv;

    // a gray pixel => set u and v zero and bail out (division by zero prevention)
    if( hsvSat<bite){
        u = v = 0.0f;
        return;
    }
    // calc classic YUV luminance y
    float y = yuv_wr*r + yuv_wg*g + yuv_wb*b;

    // calc classic YUV chrominance u and v and normalize them
    u = b-y;
    u *= hsvValInv * yuv_un;
    v = r-y;
    v *= hsvValInv * yuv_vn;

    // now we have a (u,v) vector pointing into a nice YUV hexagon in the u-v plane
    // now make the hexagon a circle
    float circleFac = hsvSat / sqrtf(u*u+v*v);  // evillive check if M4 FPU is used
    u *= circleFac;
    v *= circleFac;
}

float ExperimentalSignature::hsvValMin() const
{
    return m_hsvValMin;
}

void ExperimentalSignature::setHsvValMin(float hsvValMin)
{
    m_hsvValMin = hsvValMin;
}

float ExperimentalSignature::hsvValMax() const
{
    return m_hsvValMax;
}

void ExperimentalSignature::setHsvValMax(float hsvValMax)
{
    m_hsvValMax = hsvValMax;
}

void ExperimentalSignature::setIsActive(bool isActive)
{
    m_isActive = isActive;
}


float ExperimentalSignature::hsvHueRange() const
{
    return m_hsvHueRange;
}

void ExperimentalSignature::setHsvHueRange(float hsvHueRange)
{
    m_hsvHueRange = hsvHueRange<hueDeltaLim ? hsvHueRange*d2r : pi;
    m_hsvCosDeltaHueMin = cosf(m_hsvHueRange);
}


float ExperimentalSignature::hsvSatMin() const
{
    return m_hsvSatMin;
}
void ExperimentalSignature::setHsvSatMin(float hsvSatMin)
{
    m_hsvSatMin = hsvSatMin;
}
float ExperimentalSignature::hsvSatMax() const
{
    return m_hsvSatMax;
}
void ExperimentalSignature::setHsvSatMax(float hsvSatMax)
{
    m_hsvSatMax = hsvSatMax;
}

void ExperimentalSignature::setPosUV(const ExpSigPos &posUV)
{
    m_posUV = posUV;
    m_hsvSatMed = sqrtf(uMed()*uMed()+vMed()*vMed());
}





