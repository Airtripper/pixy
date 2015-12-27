#ifndef EXPERIMENTALSIGNATURE_H
#define EXPERIMENTALSIGNATURE_H

#include <stdint.h>
#include <cstring>
#include <math.h>

#ifdef PIXY
#include "misc.h"
#endif

class IterPixel;

// some constants
const float rgbNorm = 1.0f / 255.0f;
const float yuv_wr = 0.299f;  // classic YUV params see wikipedia
const float yuv_wb = 0.114f;
const float yuv_wg = 1.0f - yuv_wr - yuv_wb;
const float yuv_un = 1.0f/(1.0f-yuv_wb);
const float yuv_vn = 1.0f/(1.0f-yuv_wr);
const float bite = 1e-6f; // something small for float comparisons
const float pi = 3.1415926536f;
const float d2r = pi/180.0f;
const float r2d = 180.0f/pi;
const float hueDeltaLim = 23.0f; // a hue delta above x deg? => deactivate the limit and accept the full hue circle

// parameter names or their format strings
extern const char* parName_eSigUse;
extern const char* parName_eSigPos;
extern const char* parName_eSigAct;
extern const char* parName_eSigHueRng;
extern const char* parName_eSigSatMin;
extern const char* parName_eSigSatMax;
extern const char* parName_eSigValMin;
extern const char* parName_eSigValMax;

#ifdef PIXY
#define SQRTF vsqrtf
#else
#define SQRTF sqrtf
#endif



// The color signatures median (u,v) position in the chroma plane
// Used in ExperimentalSignature class below.
// Separated for storage as hidden parameter.
struct ExpSigPos{
    ExpSigPos():m_uMed(0.0f),m_vMed(0.0f){}
    float m_uMed;
    float m_vMed;
};

/// experimental signatures
/// the signature's color model related stuff
class ExperimentalSignature{

    /// median (u,v) position of the signature training set
    ExpSigPos m_posUV;
    float m_hsvSatMed; // hsv saturation of the signatures median anchor, in fact its its length |u,v|
    /// acceptance ranges (HSV)
    float m_hsvValMin;
    float m_hsvValMax;
    float m_hsvSatMin;
    float m_hsvSatMax;
    float m_hsvCosDeltaHueMin;
    float m_hsvHueRange;
    bool m_isActive;

public:

    ExperimentalSignature();
    ~ExperimentalSignature();

    /// returns true if given rgb vector is compatible with this signature
    /// The interface is a fat one as this function is also used to calculate the reference parameters
    /// u[-1..1],v[-1..1],sat[0..1] and val [0..1]
    /// r,g,b [0.0,1.0]
    /// Check if signature is activate before calling!
    /// All reference parameter could in fact be local variables of the function,
    /// but these precalculated intermediate results are still useful after calling.
    /// They are valid only if the function returns "true"!
    bool isRgbAccepted( float r, float g, float b, float& u, float& v, float& hsvVal, float& hsvSat, float& dotProd) const;

    /// same as above, but with integer r, g and b parameters [0..255]
    inline bool isRgbAccepted( uint16_t r, uint16_t g, uint16_t b, float& u, float& v, float& val, float& sat, float& dotProd) const {
        return isRgbAccepted( r*rgbNorm, g*rgbNorm, b*rgbNorm,  u, v, val, sat, dotProd);
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

    float hsvValMin() const;
    void setHsvValMin(float hsvValMin);

    inline float hsvValMax() const{return m_hsvValMax;}
    void setHsvValMax(float hsvValMax);

    inline bool isActive() const {return m_isActive;}
    void setIsActive(bool isActive);

    float hsvHueRange() const;
    void setHsvHueRange(float hsvHueRange, bool omitRangeBoost=false);

    float hsvSatMin() const;
    void setHsvSatMin(float hsvSatMin);

    float hsvSatMax() const;
    void setHsvSatMax(float hsvSatMax);

    inline float hsvSatMed() const {return m_hsvSatMed;}

    inline float uMed() const {return m_posUV.m_uMed;}
    inline float vMed() const {return m_posUV.m_vMed;}

    //ExpSigPos& accessPosUV() {return m_posUV;}
    const ExpSigPos& posUV() const{return m_posUV;}
    void setPosUV(const ExpSigPos &posUV);

};


/// A simple histogram class
class Histo{

    uint16_t* m_bins;   // the bins on heap
    uint16_t m_nBins;   // number of bins
    uint16_t m_n;       // number of entries
    uint16_t m_low;     // number of entries below bin range
    float m_min;        // lower bin range limit
    float m_max;        // upper bin range limit
    float m_binWidth;   // bin width
    float m_binWidthInv;// inverse bin width
    float m_sum;        // sum of entries
    float m_sum2;       // sum of squared etries

    /// private as cpy ctor and assignment isn't implemented
    Histo(Histo&);
    Histo& operator=(Histo&);

public:

    Histo( float min=0.0f, float max=1.0f, uint16_t nBins=32);
    ~Histo(){ delete[] m_bins; }

    /// reset histo to new range
    void reset(float min=0.0f, float max=0.0f);

    void add(float val);

    //float binVal(uint_16 b)const{ return m_min+(b+0.5f)*m_binWidth;}
    //uint16_t binCnt(uint16_t b)const{ return m_bins[b];}
    float binWidth()const{return m_binWidth;}
    float min()const{return m_min;}
    float max()const{return m_max;}
    uint16_t n()const{return m_n;}


    float mean()const{ return m_sum/m_n;}
    /// returns sigma calculated from sum and sum of squares. Might be inaccurate!
    float sigma()const;

    /// get value that corresponds to given value [0,1] of the cumulative distribution
    float X(float phi)const;
    /// get the cumulative distribution Phi for a given value X
    float phi(float X)const;
};

#endif // EXPERIMENTALSIGNATURE_H
