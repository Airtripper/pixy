#include "experimentalsignature.h"
#include "colorlut.h"  // cyclic due to IterPixel

#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;



ExperimentalSignature::ExperimentalSignature():
    m_isActive(false),
    m_uMed(0.0f),
    m_vMed(0.0f),
    m_hsvHueRange(0.0f),
    m_hsvSatRange(0.0f),
    m_hsvValMin(0.0f),
    m_hsvValMax(1.0f),
    m_hsvSatRangeFac(0.0f),
    m_cosHueRangeFac(0.0f),
    m_hsvSatMed(0.0f),
    m_hsvCosDeltaHueMin(1.0f),
    m_hsvSatMin(1.0f),
    m_hsvSatMax(0.0f)
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

    if(!m_isActive)return false;
    // calc rgb min and max
    float hsvVal = g; // aka rgbMax
    float rgbMin = g;
    if(r>hsvVal)hsvVal=r; else rgbMin=r;
    if(b>hsvVal)hsvVal=b;
    if(b<rgbMin)rgbMin=b;

    // bail out early if not within hsv value limits
    // and a hopefully correct handling for black pixels (division by zero prevention)
    if( hsvVal<m_hsvValMin+bite || hsvVal>m_hsvValMax ){
        u=v=0;
        return hsvVal+m_hsvValMin+m_hsvSatMin<bite;
    }

    float hsvValInv=1.0f/hsvVal;
    float hsvSat = (hsvVal-rgbMin) * hsvValInv;

    // bail out early if not within hsv saturation limits
    // and a hopefully correct handling for grey pixels (division by zero prevention)
    if( hsvSat<m_hsvSatMin+bite || hsvSat>m_hsvSatMax ) {
        u=v=0;
        return hsvSat+m_hsvSatMin<bite ;
    }

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
    float circleFac = hsvSat / sqrt(u*u+v*v);  // evillive check if M4 FPU is used
    u *= circleFac;
    v *= circleFac;

    // calculate cosine delta hue using dot product: (u,v)_sig . (u,v)_pix / |(u,v)_sig| / |(u,v)_pix|
    float cosHue = u*m_uMed + v*m_vMed;
    cosHue /= hsvSat * m_hsvSatMed;

    //if(cosHue>0.99) EXPLOG("u=%.2f v=%.2f um=%.2f vm=%.2f c=%.24f ", u,v,m_uMed ,m_vMed,cosHue);

    // Are we within HSV hue limits?
    return cosHue>m_hsvCosDeltaHueMin;
}

void ExperimentalSignature::init( IterPixel& pixIter)
{
    // one empty loop to get the number of pixels
    size_t sz=0;
    RGBPixel rgbPix;
    pixIter.reset();
    while(pixIter.next( 0 ,&rgbPix)) ++sz;

    // fill two vectors with all u and v values
    // Guess the M4 will have a memory problem here
    // Maybe the binary search in the ColorLUT "iterate" chain can be abused to calculate the medians
    std::vector<float> uVec(sz);
    std::vector<float> vVec(sz);
    size_t pos=0;
    pixIter.reset();
    while(pixIter.next( 0 ,&rgbPix)){
         float hsvSat, hsvVal;
         calcUV( (int16_t)rgbPix.m_r , (int16_t)rgbPix.m_g, (int16_t)rgbPix.m_b, uVec[pos], vVec[pos], hsvSat, hsvVal );
         ++pos;
    }
    // calc the u and v median, don't care if sz is odd or even
    size_t half = sz/2;
    std::nth_element( uVec.begin(), uVec.begin()+half, uVec.end());
    m_uMed = uVec[half];
    std::nth_element( vVec.begin(), vVec.begin()+half, vVec.end());
    m_vMed = vVec[half];

    // setup the precalculated helpers
    m_hsvSatMed = sqrt(m_uMed*m_uMed + m_vMed*m_vMed);

    // now investigate the saturation and hue distribution
    // to calculate reasonable hue and sat range factors
    // we are abusing the vectors used before
    pos=0;
    pixIter.reset();
    while(pixIter.next( 0 ,&rgbPix)){
         float hsvSat, hsvVal, u, v;
         calcUV( (int16_t)rgbPix.m_r , (int16_t)rgbPix.m_g, (int16_t)rgbPix.m_b, u, v, hsvSat, hsvVal );
         uVec[pos] = fabs(hsvSat-m_hsvSatMed);
         if(hsvSat>bite){
             float cosDeltaHue = u*m_uMed + v*m_vMed;
             cosDeltaHue /= m_hsvSatMed*hsvSat;
             vVec[pos] = 1.0f-cosDeltaHue;
         }else{
             vVec[pos]=0.0f;
         }
         ++pos;
    }

    std::nth_element( uVec.begin(), uVec.begin()+half, uVec.end());
    m_hsvSatRangeFac = uVec[half];
    std::nth_element( vVec.begin(), vVec.begin()+half, vVec.end());
    m_cosHueRangeFac = vVec[half];

    m_hsvSatRangeFac *= 2.0f;
    m_cosHueRangeFac *= 2.0f;

    if(m_hsvSatRangeFac<0.01f)m_hsvSatRangeFac=0.01f;
    if(m_cosHueRangeFac<0.001f)m_cosHueRangeFac=0.001f;

    // update saturation/hue filter helpers
    setHsvSatRange(m_hsvSatRange);
    setHsvHueRange(m_hsvHueRange);

    EXPLOG("init: u_median=%.2f v_median=%.2f sat_(u,v)=%.2f satRangeFac=%.2f cosHueRangeFac=%.4f", m_uMed, m_vMed, m_hsvSatMed, m_hsvSatRangeFac, m_cosHueRangeFac);
    m_isActive = true;
}

void ExperimentalSignature::calcUV(float r, float g, float b, float &u, float &v, float& hsvSat, float& hsvVal) const
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
    float circleFac = hsvSat / sqrt(u*u+v*v);  // evillive check if M4 FPU is used
    u *= circleFac;
    v *= circleFac;
}

float ExperimentalSignature::hsvSatRange() const
{
    return m_hsvSatRange;
}

void ExperimentalSignature::setHsvSatRange(float hsvSatRange)
{
    m_hsvSatRange = hsvSatRange;
    const float mx = max(m_hsvSatMed, 1.0f-m_hsvSatMed);
    const float fac = pow(mx/m_hsvSatRangeFac, m_hsvSatRange);
    EXPLOG("SatRng: %.2f %f %f",m_hsvSatRange, fac*m_hsvSatRangeFac, m_hsvSatMed);
    m_hsvSatMin = m_hsvSatMed - m_hsvSatRangeFac * fac;
    if(m_hsvSatMin<0.0f) m_hsvSatMin=0.0f;
    m_hsvSatMax = m_hsvSatMed + m_hsvSatRangeFac * fac;
    if(m_hsvSatMax>1.0f) m_hsvSatMax=1.0f;
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

bool ExperimentalSignature::isActive() const
{
    return m_isActive;
}

void ExperimentalSignature::setIsActive(bool isActive)
{
    m_isActive = isActive;
}
float ExperimentalSignature::uMed() const
{
    return m_uMed;
}

float ExperimentalSignature::vMed() const
{
    return m_vMed;
}
float ExperimentalSignature::hsvHueRange() const
{
    return m_hsvHueRange;
}

void ExperimentalSignature::setHsvHueRange(float hsvHueRange)
{
    m_hsvHueRange = hsvHueRange;
    const float fac = pow(2.0f/m_cosHueRangeFac, m_hsvHueRange);
    m_hsvCosDeltaHueMin = 1.0f - m_cosHueRangeFac * fac;
    if(m_hsvCosDeltaHueMin>1.0f) m_hsvCosDeltaHueMin=1.0f;
    if(m_hsvCosDeltaHueMin<-1.0f) m_hsvCosDeltaHueMin=-1.0f;
    EXPLOG("HueRng: %.2f %f %f %f", hsvHueRange, fac, m_cosHueRangeFac, m_hsvCosDeltaHueMin);
}



