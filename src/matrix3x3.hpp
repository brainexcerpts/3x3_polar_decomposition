#ifndef TOOL_BOX_MATRIX3X3_HPP
#define TOOL_BOX_MATRIX3X3_HPP

#include <algorithm>
#include <cassert>
#include "vec3.hpp"

/**
 * @name Mat3
 * @brief Handling 3*3 matrix
 *
 * Memory layout
  =============
  The 3x3 matrix is stored linearly with 9 floats.
  The matrix is a <b>row first layout (i.e row major)</b>

 */
struct Mat3 {

    // row major
    float a, b, c; ///< first row
    float d, e, f; ///< second row
    float g, h, i; ///< third row

    Mat3() {   }

    Mat3(const Mat3& m) {
        a = m.a; b = m.b; c = m.c;
        d = m.d; e = m.e; f = m.f;
        g = m.g; h = m.h; i = m.i;
    }

    Mat3(float a_, float b_, float c_,
         float d_, float e_, float f_,
         float g_, float h_, float i_)
    {
        a = a_; b = b_; c = c_;
        d = d_; e = e_; f = f_;
        g = g_; h = h_; i = i_;
    }

    Mat3(const Vec3& x,
         const Vec3& y,
         const Vec3& z)
    {
        a = x.x; b = y.x; c = z.x;
        d = x.y; e = y.y; f = z.y;
        g = x.z; h = y.z; i = z.z;
    }

    static Mat3 diagonal(float d)
    {

        return Mat3(d   , 0.0f, 0.0f,
                    0.0f, d   , 0.0f,
                    0.0f, 0.0f, d );

    }

    static inline Mat3 from_array_ptr( const float* row_major_ptr){
        Mat3 result;
        for(int i = 0; i < 9; ++i)
            result[i] = row_major_ptr[i];
        return result;
    }

    // -------------------------------------------------------------------------
    /// @name operators overload
    // -------------------------------------------------------------------------


    //----------------
    // Multiplications
    //----------------

    Vec3 operator*(const Vec3& v) const
    {
        float x = v.x * a + v.y * b + v.z * c;
        float y = v.x * d + v.y * e + v.z * f;
        float z = v.x * g + v.y * h + v.z * i;
        return Vec3(x, y, z);
    }

    inline Mat3 operator*(const Mat3& m) const
    {
        return Mat3(a * m.a + b * m.d + c * m.g,
                         a * m.b + b * m.e + c * m.h,
                         a * m.c + b * m.f + c * m.i,
                         d * m.a + e * m.d + f * m.g,
                         d * m.b + e * m.e + f * m.h,
                         d * m.c + e * m.f + f * m.i,
                         g * m.a + h * m.d + i * m.g,
                         g * m.b + h * m.e + i * m.h,
                         g * m.c + h * m.f + i * m.i);
    }

    inline Mat3 operator*(float x) const
    {
        return Mat3(a * x, b * x, c * x,
                         d * x, e * x, f * x,
                         g * x, h * x, i * x);
    }

    Mat3& operator*=(float x)
    {
        a *= x; b *= x; c *= x;
        d *= x; e *= x; f *= x;
        g *= x; h *= x; i *= x;
        return *this;
    }

    friend
    Mat3 operator*(const float x_, const Mat3& mat)
    {
        return Mat3(x_ * mat.a, x_ * mat.b, x_ * mat.c,
                         x_ * mat.d, x_ * mat.e, x_ * mat.f,
                         x_ * mat.g, x_ * mat.h, x_ * mat.i);
    }

    //----------
    // Divisions
    //----------

    Mat3 operator/(float x) const
    {
        return Mat3(a / x, b / x, c / x,
                         d / x, e / x, f / x,
                         g / x, h / x, i / x);
    }

private:
    // Dividing a matrix by another doesn't have any clear meaning
    Mat3 operator/(const Mat3& m) /* = delete*/ const;
public:

    Mat3& operator/=(float x)
    {
        a /= x; b /= x; c /= x;
        d /= x; e /= x; f /= x;
        g /= x; h /= x; i /= x;
        return *this;
    }

    friend Mat3 operator/(const float x_, const Mat3& mat)
    {
        return Mat3(x_ / mat.a, x_ / mat.b, x_ / mat.c,
                         x_ / mat.d, x_ / mat.e, x_ / mat.f,
                         x_ / mat.g, x_ / mat.h, x_ / mat.i);
    }

    //----------
    // Additions
    //----------


    Mat3 operator+(const Mat3& m) const
    {
        return Mat3(a + m.a, b + m.b, c + m.c,
                         d + m.d, e + m.e, f + m.f,
                         g + m.g, h + m.h, i + m.i);
    }

    Mat3 operator+(float x) const
    {
        return Mat3(a + x, b + x, c + x,
                         d + x, e + x, f + x,
                         g + x, h + x, i + x);
    }

    friend Mat3 operator+(const float x_, const Mat3& mat)
    {
        return Mat3(x_ + mat.a, x_ + mat.b, x_ + mat.c,
                         x_ + mat.d, x_ + mat.e, x_ + mat.f,
                         x_ + mat.g, x_ + mat.h, x_ + mat.i);
    }

    Mat3& operator+=(float x)
    {
        a += x; b += x; c += x;
        d += x; e += x; f += x;
        g += x; h += x; i += x;
        return *this;
    }

    //--------------
    // Substractions
    //--------------

    Mat3 operator-(const Mat3& m) const
    {
        return Mat3(a - m.a, b - m.b, c - m.c,
                         d - m.d, e - m.e, f - m.f,
                         g - m.g, h - m.h, i - m.i);
    }

    Mat3 operator-() const
    {
        return Mat3(-a, -b, -c,
                         -d, -e, -f,
                         -g, -h, -i);
    }

    Mat3 operator-(float x) const
    {
        return Mat3(a - x, b - x, c - x,
                         d - x, e - x, f - x,
                         g - x, h - x, i - x);
    }

    friend Mat3 operator-(const float x_, const Mat3& mat)
    {
        return Mat3(x_ - mat.a, x_ - mat.b, x_ - mat.c,
                         x_ - mat.d, x_ - mat.e, x_ - mat.f,
                         x_ - mat.g, x_ - mat.h, x_ - mat.i);
    }

    Mat3& operator-=(float x)
    {
        a -= x; b -= x; c -= x;
        d -= x; e -= x; f -= x;
        g -= x; h -= x; i -= x;
        return *this;
    }

    // -------------------------------------------------------------------------
    /// @name Accessors
    // -------------------------------------------------------------------------

    /// Conversion returns the memory address of the matrix.
    /// (row major)
    explicit operator const float*() const { return ((float*)(this)); }
    /// Conversion returns the memory address of the vector. (Non const version)
    explicit operator float*() { return ((float*)(this)); }

    /// Conversion returns the memory address of the matrix.
    /// (row major)
    const float* data() const { return ((float*)(this)); }
    /// Conversion returns the memory address of the vector. (Non const version)
    float* data() { return ((float*)(this)); }

    //----------------
    // Access elements
    //----------------

    inline
    const float& operator()(int row, int column) const
    {
        assert(row >= 0 && row < 3);
        assert(column >= 0 && column < 3);
        return ((float*)(this))[column + row*3];
    }

    inline
    float& operator()(int row, int column)
    {
        assert(row >= 0 && row < 3);
        assert(column >= 0 && column < 3);
        return ((float*)(this))[column + row*3];
    }

    float& operator[](int idx) {
        assert(idx >= 0 && idx < 9);
        return ((float*)(this))[idx];
    }

    const float& operator[](int idx) const {
        assert(idx >= 0 && idx < 9);
        return ((float*)(this))[idx];
    }

    // -------------------------------------------------------------------------
    /// @name operations
    // -------------------------------------------------------------------------

    float det() const
    {
        return a * ( e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    }

    /// @return the matrix with normalized x, y, z column vectors
    /// (basically eliminates scale factors of the matrix)
    Mat3 normalized() const {
        return Mat3(x().normalized(), y().normalized(), z().normalized());
    }

    Mat3 inverse() const
    {
        float c0 = e * i - f * h;
        float c1 = f * g - d * i;
        float c2 = d * h - e * g;
        float idet = (a * c0 + b * c1 + c * c2);
        if(false){
            idet = float(1.0) / idet;
            return Mat3(c0 , c * h - b * i, b * f - c * e,
                             c1 , a * i - c * g, c * d - a * f,
                             c2 , b * g - a * h, a * e - b * d) * idet;
        }else{
            return Mat3(c0 , c * h - b * i, b * f - c * e,
                             c1 , a * i - c * g, c * d - a * f,
                             c2 , b * g - a * h, a * e - b * d) / idet;
        }

    }

    Mat3 transpose() const
    {
        return Mat3(a, d, g,
                         b, e, h,
                         c, f, i);
    }

    void set_abs()
    {
        a = std::abs(a); b = std::abs(b); c = std::abs(c);
        d = std::abs(d); e = std::abs(e); f = std::abs(f);
        g = std::abs(g); h = std::abs(h); i = std::abs(i);
    }


    float max_elt() const
    {
        return std::max(i, std::max(std::max(std::max(a,b),std::max(c,d)),
                                    std::max(std::max(e,f),std::max(g,h))));
    }

    float min_elt() const
    {
        return std::min(i, std::min(std::min(std::min(a,b),std::min(c,d)),
                                    std::min(std::min(e,f),std::min(g,h))));
    }

    Mat3 get_ortho() const
    {
        Mat3 h0 = (*this);
        Mat3 h1 = h0;
        h1.set_abs();
        float eps = (float(1) +  h1.min_elt()) * float(1e-5);
        for(int i = 0; i < 500/* to avoid infinite loop */; i++){
            h0 = (h0 + (h0.inverse()).transpose()) * float(0.5);
            h1 = h1 - h0;
            h1.set_abs();
            if(h1.max_elt() <= eps)
                break;
            h1 = h0;
        }
        return h0;
    }

    float get_rotation_axis_angle(Vec3& axis) const
    {
        axis.x = h - f + float(1e-5);
        axis.y = c - g;
        axis.z = d - b;
        float sin_angle = axis.safe_normalize();
        float cos_angle = a + e + i - float(1);
        return std::atan2(sin_angle, cos_angle);
    }

    Vec3 x() const { return Vec3(a, d, g); }
    Vec3 y() const { return Vec3(b, e, h); }
    Vec3 z() const { return Vec3(c, f, i); }

    void set_x(const Vec3& x) { a = x.x; d = x.y; g = x.z; }
    void set_y(const Vec3& y) { b = y.x; e = y.y; h = y.z; }
    void set_z(const Vec3& z) { c = z.x; f = z.y; i = z.z; }

    inline
    void set(float a_, float b_, float c_,
             float d_, float e_, float f_,
             float g_, float h_, float i_)
    {
        a = a_; b = b_; c = c_;
        d = d_; e = e_; f = f_;
        g = g_; h = h_; i = i_;
    }

    //--------------------------------------------------------------------------
    /// @name Static constructors
    //--------------------------------------------------------------------------

    static Mat3 identity()
    {
        return Mat3(float(1), float(0), float(0),
                         float(0), float(1), float(0),
                         float(0), float(0), float(1));

    }

    /// @return the rotation matrix given 'axis' and 'angle' in radian
    static Mat3 rotate(const Vec3& axis, float radian_angle)
    {
        Vec3 n = axis;
        n.normalize();
        float cp = std::cos(radian_angle);
        float sp = std::sin(radian_angle);
        float acp = float(1) - cp;
        float nxz = n.x * n.z;
        float nxy = n.x * n.y;
        float nyz = n.y * n.z;
        return Mat3(cp + acp * n.x * n.x,
                         acp * nxy - sp * n.z,
                         acp * nxz + sp * n.y,

                         acp * nxy + sp * n.z,
                         cp + acp * n.y * n.y,
                         acp * nyz - sp * n.x,

                         acp * nxz - sp * n.y,
                         acp * nyz + sp * n.x,
                         cp + acp * n.z * n.z);
    }

};


#endif // TOOL_BOX_MATRIX3X3_HPP
