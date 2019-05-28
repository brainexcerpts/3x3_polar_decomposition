#ifndef TOOL_BOX_POLAR_DECOMPOSITION_3x3_HPP
#define TOOL_BOX_POLAR_DECOMPOSITION_3x3_HPP

#include "matrix3x3.hpp"


/// @note : code adapted from Vega 3.0 (http://run.usc.edu/vega/download.html)
/// @tparam force_rotation :
/// if true, we further ensure the matrix R will be a rotation,
/// S will be symmetric, but not necessarily positive-definite.
/// (this triggers an extra computation of a 3x3 SVD when the standard
/// routine is unable to perform the decomposition
template<bool force_rotation>
class Polar_decomposition {
public:

    Polar_decomposition() {}

    /// @param[in] M : is 3x3 input matrix (linear array [9], row major,
    /// tbx::Mat3 is also a valid argument)
    Polar_decomposition(const Mat3& M, float tolerance = 1E-6){
        compute(M, tolerance);
    }

    /// Computes the Polar Decomposition of a general 3x3 matrix M.
    /// M = R * S
    /// @param[in] M : is 3x3 input matrix (linear array [9], row major,
    /// tbx::Mat3 is also a valid argument)
    ///
    /// Compute will caculate:
    /// R is 3x3 orthogonal matrix, R R^T = R^T R = I (rotation)
    /// S is 3x3 symmetric positive-definite matrix
    /// @return det(R) = sgn(det(M)); this sign can be 1 or -1, depending on M
    float compute(const Mat3& M, float tolerance = 1E-6);

    const Mat3& matrix_R() const { return _R; }
    const Mat3& matrix_S() const { return _S; }

private:
    Mat3 _R; ///< Rotation
    Mat3 _S; ///< Scale and shear
};

#include "polar_decomposition_3x3.inl"

#endif // TOOL_BOX_POLAR_DECOMPOSITION_3x3_HPP
