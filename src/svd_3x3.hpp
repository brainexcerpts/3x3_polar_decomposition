#ifndef TOOL_BOX_SVD_3x3_HPP
#define TOOL_BOX_SVD_3x3_HPP

#include "matrix3x3.hpp"

// Note : code adapted from Vega 3.0 http://run.usc.edu/vega/download.html
// see: libraries/minivector/mat3d.h -> void SVD(...)
// (released under the 3-clause BSD license,
// which means that it can be used freely both in academic
// research and in commercial applications. )

/**
  Standard SVD (modifiedSVD == 0):

  Given a 3x3 matrix M, decomposes it using SVD so
  that M = U Sigma V^T where U is a 3x3 rotation,
  V is 3x3 orthonormal (V^T V = V V^T = I), and
  Sigma is a diagonal matrix with non-negative entries,
  in descending order, Sigma[0] >= Sigma[1] >= Sigma[2] >= 0.
  Note that the matrix V may have a determinant equaling
  to -1 (i.e., reflection).

  singularValue_eps is a threshold to determine
  when a singular value is deemed zero, and special handling is then invoked
  to improve numerical robustness.

  Modified SVD (modifiedSVD == 1):

  The SVD is modified so that it gives the
  following properties (useful in solid mechanics applications) :

  1) Not just the determinant of U, but also the determinant of V is 1
    (i.e., both U and V are rotations, not reflections).

  2) The smallest diagonal value Sigma[2] may be negative. We have:
     Sigma[0] >= Sigma[1] >= abs(Sigma[2]) >= 0 .

    ----------------------------------------------------------

    The implementation follows section 5 of
    G. Irving, J. Teran, and R. Fedkiw.
    Invertible Finite Elements for Robust Simulation of Large Deformation.
    In Symp. on Computer Animation (SCA), pages 131-140, 2004.
    It computes a^T a, and computes its eigenvectors, a^T a = V Sigma^2 V^T.
    Sigma is then recovered using sqrt from Sigma^2.
    To recover U, compute U = F * V * diag(Sigma^{-1}).
    Care must be taken when singular values of Sigma are small (this is
    handled in the code).

    We have also added the ability to compute standard SVD (modifiedSVD = 0).
*/
int SVD_3x3(const Mat3& F, Mat3& U, Vec3& sigma, Mat3& V, float eps, int modifiedSVD);

#include "svd_3x3.inl"

#endif // TOOL_BOX_SVD_3x3_HPP
