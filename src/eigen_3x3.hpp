#ifndef TOOL_BOX_EIGEN_3X3_HPP
#define TOOL_BOX_EIGEN_3X3_HPP

/// @brief Eigen-decomposition for symmetric 3x3 real matrices.
/// @note Public domain

/// Computes eigenvalues and eigenvectors of a 3x3 matrix M
/// Assumes symmetric matrix;
/// contents of matrix "M" are not modified by the routine
/// Eigenvalues are sorted in decreasing order
/// (not decreasing absolute-value order)
/// Returned eigenvectors are unit length
/// Symmetric matrix A => eigenvectors in columns of V, corresponding
/// eigenvalues in d.
void eigen_decomposition(const float A[3][3], float V[3][3], float d[3]);

#endif // TOOL_BOX_EIGEN_3X3_HPP
