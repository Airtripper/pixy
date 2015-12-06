#ifndef EXPERIMENTALSIGNATURE_H
#define EXPERIMENTALSIGNATURE_H

#include <stdint.h>

class IterPixel;

// some constants
const float rgbNorm = 1.0f / 255.0f;
const float yuv_wr = 0.299f;
const float yuv_wb = 0.114f;
const float yuv_wg = 1.0f - yuv_wr - yuv_wb;
const float yuv_un = 1.0f/(1.0f-yuv_wb);
const float yuv_vn = 1.0f/(1.0f-yuv_wr);
// something small for float comparisons
const float bite = 1e-6f;

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
    void calcUV(float r, float g, float b, float& u, float& v, float& hsvSat, float& hsvVal) const;

    /// same as above, but with integer r, g and b parameters [0..255]
    inline void calcUV( int16_t r, int16_t g, int16_t b, float& u, float& v, float& hsvSat, float& hsvVal) const {
        return calcUV( r*rgbNorm, g*rgbNorm, b*rgbNorm,  u, v, hsvSat, hsvVal);
    }


public:
    float hsvSatRange() const;
    void setHsvSatRange(float hsvSatRange);

    float hsvValMin() const;
    void setHsvValMin(float hsvValMin);

    float hsvValMax() const;
    void setHsvValMax(float hsvValMax);

    bool isActive() const;
    void setIsActive(bool isActive);

    float uMed() const;
    float vMed() const;

    float hsvHueRange() const;
    void setHsvHueRange(float hsvHueRange);

private:
    bool m_isActive;
    // median position of the signature training set
    float m_uMed;
    float m_vMed;
    // acceptance ranges (HSV)
    float m_hsvHueRange;
    float m_hsvSatRange;
    float m_hsvValMin;
    float m_hsvValMax;
    // range factors (slider value translation)
    float m_hsvSatRangeFac;
    float m_cosHueRangeFac;
    // precalculated helpers
    float m_hsvSatMed; // hsv saturation of the signatures median anchor, in fact its its length |u,v|
    float m_hsvCosDeltaHueMin;
    float m_hsvSatMin;
    float m_hsvSatMax;
};

#endif // EXPERIMENTALSIGNATURE_H
