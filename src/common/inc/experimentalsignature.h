#ifndef EXPERIMENTALSIGNATURE_H
#define EXPERIMENTALSIGNATURE_H

#include <stdint.h>
#include <cstring>
#include <math.h>

class IterPixel;

// some constants
const float rgbNorm = 1.0f / 255.0f;
const float yuv_wr = 0.299f;
const float yuv_wb = 0.114f;
const float yuv_wg = 1.0f - yuv_wr - yuv_wb;
const float yuv_un = 1.0f/(1.0f-yuv_wb);
const float yuv_vn = 1.0f/(1.0f-yuv_wr);
const float bite = 1e-6f; // something small for float comparisons
const float pi = 3.1415926536f;
const float d2r = pi/180.0f;
const float r2d = 180.0f/pi;

// The color signatures median (u,v) position in the chroma plane
// Used in ExperimentalSignature class below.
// Separated for storage as hidden parameter.
struct ExpSigPos{
    float m_uMed;
    float m_vMed;
};

/// experimental signatures
/// the signature's color model related stuff
class ExperimentalSignature
{
public:

    ExperimentalSignature();
    ~ExperimentalSignature();

    /// returns true if given rgb vector is compatible with this signature
    /// The interface is a fat one as this function is also used to calculate the reference parameters
    /// u[-1..1],v[-1..1],sat[0..1] and val [0..1]
    /// r,g,b [0.0,1.0]
    bool isRgbAccepted( float r, float g, float b, float& u, float& v) const;

    /// same as above, but with integer r, g and b parameters [0..255]
    inline bool isRgbAccepted( int16_t r, int16_t g, int16_t b, float& u, float& v) const {
        return m_isActive && isRgbAccepted( r*rgbNorm, g*rgbNorm, b*rgbNorm,  u, v);
    }

    /// initialize this signature from the given pixel iterator
    void init( IterPixel& pixIter );

    /// a helper that calculates the u and v chroma valus from the given rgb vector
    /// u[-1..1],v[-1..1],sat[0..1] and val [0..1]
    /// r,g,b [0.0,1.0]
    void static translateRGB(float r, float g, float b, float& u, float& v, float& hsvSat, float& hsvVal) ;

    /// same as above, but with integer r, g and b parameters [0..255]
    inline static void translateRGB( uint8_t r, uint8_t g, uint8_t b, float& u, float& v, float& hsvSat, float& hsvVal)  {
        return translateRGB( r*rgbNorm, g*rgbNorm, b*rgbNorm,  u, v, hsvSat, hsvVal);
    }


public:

    float hsvValMin() const;
    void setHsvValMin(float hsvValMin);

    float hsvValMax() const;
    void setHsvValMax(float hsvValMax);

    bool isActive() const;
    void setIsActive(bool isActive);

    float hsvHueRange() const;
    void setHsvHueRange(float hsvHueRange);

    float hsvSatMin() const;
    void setHsvSatMin(float hsvSatMin);

    float hsvSatMax() const;
    void setHsvSatMax(float hsvSatMax);

    inline float uMed() const {return m_posUV.m_uMed;}
    inline float vMed() const {return m_posUV.m_vMed;}

    //ExpSigPos& accessPosUV() {return m_posUV;}
    const ExpSigPos& posUV() const{return m_posUV;}
    void setPosUV(const ExpSigPos &posUV);

private:
    // median (u,v) position of the signature training set
    ExpSigPos m_posUV;
    float m_hsvSatMed; // hsv saturation of the signatures median anchor, in fact its its length |u,v|
    // acceptance ranges (HSV)
    float m_hsvValMin;
    float m_hsvValMax;
    float m_hsvSatMin;
    float m_hsvSatMax;
    float m_hsvCosDeltaHueMin;
    float m_hsvHueRange;
    bool m_isActive;
};


/// A simple histogram class
struct Histo{

    static const uint16_t s_nBins=32;
    uint16_t m_bins[s_nBins];
    uint16_t m_low;
    float m_min;
    float m_max;
    float m_binWidth;
    float m_binWidthInv;
    float m_sum;
    float m_sum2;

    uint16_t m_n;

    void reset(float min, float max){
        memset( m_bins, 0, sizeof(uint16_t)*s_nBins);
        m_low = 0;
        m_min = min;
        m_max = max;
        m_binWidth = (max-min)/s_nBins;
        m_binWidthInv = 1.0f/m_binWidth;
        m_sum = m_sum2 = 0;
        m_n=0;
    }

    Histo( float min, float max){ reset(min,max);}

    void add(float val){
        ++m_n;
        m_sum+=val;
        m_sum2+=val*val;
        if(val<m_min) ++m_low;
        else{
            uint16_t bin = (val-m_min)*m_binWidthInv;
            if(bin<s_nBins) ++m_bins[bin];
        }
    }

    // float binVal(uint_16 b){ return m_min+(b+0.5f)*m_binWidth;}
    // uint16_t binCnt(uint16_t b){ return m_bins[b];}

    float mean(){ return m_sum/m_n;}

    float sigma(){
        float var = m_sum2 - m_sum*m_sum/m_n; // acuracy? I should have read D.K.
        var /= m_n-1;
        return var>0 ? sqrt(var) : 0.0f;
    }

    /// get value that corresponds to given value [0,1] of the cumulative distribution
    float X(float phi){
        uint16_t lim=phi*m_n+0.5f;
        uint16_t sum = m_low;
        if(sum>lim) return m_min;
        for(uint16_t b=0; b<s_nBins; ++b){
            if(sum+m_bins[b]>lim) return m_min + m_binWidth * ( b + float(lim-sum) / m_bins[b] );
            sum += m_bins[b];
        }
        return m_max;
    }
};

#endif // EXPERIMENTALSIGNATURE_H
