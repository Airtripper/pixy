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

    // calc rgb min and max
    float hsvVal = g; // aka rgbMax
    float rgbMin = g;
    if(r>hsvVal)hsvVal=r; else rgbMin=r;
    if(b>hsvVal)hsvVal=b;
    if(b<rgbMin)rgbMin=b;

    // bail out early if not within hsv value limits
    // and handling of black pixels (division by zero prevention)
    // The additional cut "hsvVal<=1.0f" handles hsv values above 1.0
    // This might happen due to floating point inaccuracies or the g1/g2 bayer mess.
    if( hsvVal<m_hsvValMin+bite || (hsvVal>m_hsvValMax && hsvVal<=1.0f) ){
        u=v=0.0f;
        return hsvVal+m_hsvValMin+m_hsvSatMin < bite;  // return true if signature accepts black
    }

    float hsvSat = hsvVal-rgbMin;  // ... not really as sat=(val-min)/val but I wanted to defer the division after the next bail out check

    // bail out early if not within hsv saturation limits
    // and handling of grey pixels (division by zero prevention)
    if( hsvSat<m_hsvSatMin*hsvVal+bite || hsvSat>m_hsvSatMax*hsvVal) {
        u=v=0;
        return hsvSat+m_hsvSatMin<bite ; // return true if signature accepts 255 shades of grey ... :D
    }

    float hsvValInv=1.0f/hsvVal;
    hsvSat *= hsvValInv;     // now its really a saturation value

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
#ifdef PIXY
    float circleFac = hsvSat / vsqrtf(u*u+v*v);  // evillive check if M4 FPU is used ... now it is!
#else
    float circleFac = hsvSat / sqrtf(u*u+v*v);
#endif
    u *= circleFac;
    v *= circleFac;

    // We got here with with zero signature saturation => Accept the candidate
    // This bail out prevents a division by zero below.
    // We could have bailed out earlier, but we wanna return valid u and v values
    if(m_hsvSatMed<bite)return true;

    // calculate cosine delta hue using dot product: (u,v)_sig . (u,v)_pix / |(u,v)_sig| / |(u,v)_pix|
    float cosHue = u*uMed() + v*vMed(); // missing the division? Made it a multiplication in the last check below

    // Are we within HSV hue limits?
    // return cosHue>m_hsvCosDeltaHueMin; // save one division
    return cosHue > m_hsvCosDeltaHueMin * hsvSat * m_hsvSatMed;
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
    Histo valHist( 0.0f, 1.0f);
    RGBPixel rgbPix;
    pixIter.reset();
    while(pixIter.next( 0 ,&rgbPix)){
         float u, v, sat, val;
         translateRGB( rgbPix.m_r,rgbPix.m_g,rgbPix.m_b,  u, v, sat, val );
         hueHist.add( u);
         satHist.add( v);
         valHist.add( val);
    }
    float uMean = hueHist.mean();
    float vMean = satHist.mean();

    // calculate the polar angle of the mean (u,v) vector
    float hueOff = fabs(uMean)>bite && fabs(vMean)>bite ? atan2(vMean,uMean) : 0.0f;

    // now analyse the saturation and hue distribution
    // hue is rotated by hueOff to the center of the histo before
    // ... If the atan2 stuff below is too expensive for the M4,
    // I guess, that it's just sufficient to take the median u and v
    // from above directly to determine the hue median and range.
    // Maybe there is no difference at all.
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

    // saturation and hue selection cuts are defined by the values
    // for 5% and 95% of the cumulative distribution
    const float outLim = 0.05f;
    const float outFac = 3.0f;
    float low = satHist.X( outLim);
    float hgh = satHist.X( 1.0f-outLim);
    float add = (hgh-low)*outLim*outFac;
    m_hsvSatMin = low-add >0.0f ? low-add : 0.0f;
    m_hsvSatMax = hgh+add <1.0f ? hgh+add : 1.0f;

    m_hsvHueRange = hueHist.X(1.0f-outLim);
    m_hsvHueRange += m_hsvHueRange * outLim * outFac;
    const float dMin = 3.0f*d2r; // don't set min hue lim smaller than LUT resolution
    if(m_hsvHueRange<dMin)m_hsvHueRange=dMin;
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

    // set the value limits depending on the saturation result
    if(m_hsvSatMed>0.2f){
        // we have some colour => set the value acceptance range to wide defaults
        m_hsvValMin = 0.2f;
        m_hsvValMax = 1.0f;
        m_hsvSatMax = 1.0f; // and also maximise the upper saturation cut
    }else{
        // rather grey => use value range of the training set
        low = valHist.X( outLim);
        hgh = valHist.X( 1.0f-outLim);
        add = (hgh-low)*outLim * outFac;
        m_hsvValMin = low-add >0.0f ? low-add : 0.0f;
        m_hsvValMax = hgh+add <1.0f ? hgh+add : 1.0f;
        m_hsvSatMin = 0.0f; // and also minimise the lower saturation cut
        if(m_hsvSatMax<0.05f)m_hsvSatMax=0.05f; // and set sat max above the LUT resolution
        setHsvHueRange(hueDeltaLim); // and set hue lim to full circle
    }

    // and calculate the u, v and hue medians
    float hueMed=hueHist.X(0.5)+hueOff;
    m_posUV.m_uMed = m_hsvSatMed*cosf( hueMed);
    m_posUV.m_vMed = m_hsvSatMed*sinf( hueMed);

    /*DBG("uMed=%.2f vMed=%.2f satMed=%.2f satMin=%.2f satMax=%.2f hueMed=%.2f hueRng=%.1f hueCos=%.4f hueOff=%.2f",
           uMed(), vMed(), m_hsvSatMed, m_hsvSatMin, m_hsvSatMax, hueMed*r2d, m_hsvHueRange*r2d, m_hsvCosDeltaHueMin, hueOff*r2d );*/

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
#ifdef PIXY
    float circleFac = hsvSat / vsqrtf(u*u+v*v);  // evillive check if M4 FPU is used ... now it is!
#else
    float circleFac = hsvSat / sqrtf(u*u+v*v);
#endif
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
    return m_hsvHueRange*r2d;
}

void ExperimentalSignature::setHsvHueRange(float hsvHueRange)
{
    m_hsvHueRange = hsvHueRange<=hueDeltaLim-bite ? hsvHueRange*d2r : pi;
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
#ifdef PIXY
    m_hsvSatMed = vsqrtf(uMed()*uMed()+vMed()*vMed());
#else
    m_hsvSatMed = sqrtf(uMed()*uMed()+vMed()*vMed());
#endif

}





void Histo::reset(float min, float max){
    memset( m_bins, 0, sizeof(uint16_t)*m_nBins);
    m_low = m_sum = m_sum2 = m_n=0;
    if(max!=min){
        m_min = min;
        m_max = max;
        m_binWidth = (max-min)/m_nBins;
        m_binWidthInv = 1.0f/m_binWidth;
    }

}

Histo::Histo( float min, float max, uint16_t nBins){
    m_nBins=nBins;
    m_bins = new uint16_t[m_nBins];
    reset( min, max);
}

void Histo::add(float val){
    ++m_n;
    m_sum+=val;
    m_sum2+=val*val;
    if(val<m_min) ++m_low;
    else{
        uint16_t bin = (val-m_min)*m_binWidthInv;
        if(bin<m_nBins) ++m_bins[bin];
    }
}

float Histo::sigma() const {
    float var = m_sum2 - m_sum*m_sum/m_n; // acuracy? I should have read D.K.
    var /= m_n-1;
#ifdef PIXY
    return var>0 ? vsqrtf(var) : 0.0f;
#else
    return var>0 ? sqrtf(var) : 0.0f;
#endif
}

float Histo::X(float phi)const{
    uint16_t lim=phi*m_n+0.5f;
    uint16_t sum = m_low;
    if(sum>lim) return m_min;
    for(uint16_t b=0; b<m_nBins; ++b){
        if(sum+m_bins[b]>lim) return m_min + m_binWidth * ( b + float(lim-sum) / m_bins[b] );
        sum += m_bins[b];
    }
    return m_max;
}




