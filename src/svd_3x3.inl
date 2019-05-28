#include "svd_3x3.hpp"

#include <cmath>

#include "eigen_3x3.hpp"
#include "vec3.hpp"

// =============================================================================
namespace details {
// =============================================================================

/// Given an input vector v, find a unit vector that is orthogonal to it
static inline
Vec3 find_orthonormal_vec(const Vec3& vec)
{
    // find smallest abs component of v
    int smallest_idx = 0;
    for(int dim = 1; dim < 3; dim++)
        if (fabs(vec[dim]) < fabs(vec[smallest_idx]))
            smallest_idx = dim;

    Vec3 axis(0.0, 0.0, 0.0);
    axis[smallest_idx] = 1.0;

    // this cross-product will be non-zero (as long as v is not zero)
    return (vec.cross(axis)).normalized();
}

// -----------------------------------------------------------------------------

static inline
Mat3 multiply_diag_right(const Mat3& mat, Vec3& v)
{
    Mat3 result = mat;

    result(0,0) *= v[0];
    result(1,0) *= v[0];
    result(2,0) *= v[0];

    result(0,1) *= v[1];
    result(1,1) *= v[1];
    result(2,1) *= v[1];

    result(0,2) *= v[2];
    result(1,2) *= v[2];
    result(2,2) *= v[2];

    return result;
}

} // END details namespace -----------------------------------------------------

static inline
void eigen_sym(Mat3& a, Vec3& eig_val, Vec3 eig_vec[3])
{
    float A[3][3] = { {a(0,0), a(0,1), a(0,2)},
                      {a(1,0), a(1,1), a(1,2)},
                      {a(2,0), a(2,1), a(2,2)} };
    float V[3][3];
    float d[3];
    eigen_decomposition(A, V, d);

    eig_val = Vec3(d[2],d[1],d[0]);
    if(eig_vec)
    {
        eig_vec[0] = Vec3(V[0][2], V[1][2], V[2][2]);
        eig_vec[1] = Vec3(V[0][1], V[1][1], V[2][1]);
        eig_vec[2] = Vec3(V[0][0], V[1][0], V[2][0]);
    }
}

// -----------------------------------------------------------------------------

inline int SVD_3x3(const Mat3 &F, Mat3& U, Vec3& sigma, Mat3& V, float eps, int modifiedSVD)
{
    // The code handles the following special situations:

    //---------------------------------------------------------
    // 1. det(V) == -1
    //    - multiply the first column of V by -1
    //---------------------------------------------------------
    // 2. An entry of Sigma is near zero
    //---------------------------------------------------------
    // (if modifiedSVD == 1) :
    // 3. negative determinant (Tet is inverted in solid mechanics).
    //    - check if det(U) == -1
    //    - If yes, then negate the minimal element of Sigma
    //      and the corresponding column of U
    //---------------------------------------------------------

    // form F^T F and do eigendecomposition
    Mat3 normalEq = F.transpose() * F;
    Vec3 eigen_vals;
    Vec3 eigen_vecs[3];

    eigen_sym(normalEq, eigen_vals, eigen_vecs);

    V.set(eigen_vecs[0][0], eigen_vecs[1][0], eigen_vecs[2][0],
          eigen_vecs[0][1], eigen_vecs[1][1], eigen_vecs[2][1],
          eigen_vecs[0][2], eigen_vecs[1][2], eigen_vecs[2][2]);

    // Handle situation:
    // 1. det(V) == -1
    //    - multiply the first column of V by -1
    if (V.det() < 0.0)
    {
        // convert V into a rotation (multiply column 1 by -1)
        V(0,0) *= -1.0;
        V(1,0) *= -1.0;
        V(2,0) *= -1.0;
    }

    sigma[0] = (eigen_vals[0] > 0.0f) ? sqrtf(eigen_vals[0]) : 0.0f;
    sigma[1] = (eigen_vals[1] > 0.0f) ? sqrtf(eigen_vals[1]) : 0.0f;
    sigma[2] = (eigen_vals[2] > 0.0f) ? sqrtf(eigen_vals[2]) : 0.0f;

    //printf("--- Sigma ---\n");
    //printf("%G %G %G\n", Sigma[0][0], Sigma[1][1], Sigma[2][2]);

    // compute inverse of singular values
    // also check if singular values are close to zero
    Vec3 sigma_inverse;
    sigma_inverse[0] = (sigma[0] > eps) ? (1.0f / sigma[0]) : 0.0f;
    sigma_inverse[1] = (sigma[1] > eps) ? (1.0f / sigma[1]) : 0.0f;
    sigma_inverse[2] = (sigma[2] > eps) ? (1.0f / sigma[2]) : 0.0f;

    // compute U using the formula:
    // U = F * V * diag(SigmaInverse)
    U = F * V;
    U = details::multiply_diag_right(U, sigma_inverse);

    // In theory, U is now orthonormal, U^T U = U U^T = I ..
    // it may be a rotation or a reflection, depending on F.
    // But in practice, if singular values are small or zero,
    // it may not be orthonormal, so we need to fix it.
    // Handle situation:
    // 2. An entry of Sigma is near zero
    // ---------------------------------------------------------

    if ((sigma[0] < eps) && (sigma[1] < eps) && (sigma[2] < eps))
    {
        // extreme case, all singular values are small,
        // material has collapsed almost to a point
        // see [Irving 04], p. 4
        U.set(1.0, 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, 1.0);
    }
    else
    {
        // handle the case where two singular values are small,
        // but the third one is not handle it by computing
        // two (arbitrary) vectors orthogonal to the eigenvector
        // for the large singular value
        int done = 0;
        for(int dim = 0; dim < 3; dim++)
        {
            int dim_a = dim;
            int dim_b = (dim + 1) % 3;
            int dim_c = (dim + 2) % 3;
            if ((sigma[dim_b] < eps) && (sigma[dim_c] < eps))
            {
                // only the column dimA can be trusted,
                // columns dimB and dimC correspond to tiny singular values
                Vec3 vec1(U(0,dim_a), U(1,dim_a), U(2,dim_a)); // column dimA
                Vec3 vec2;
                vec2 = details::find_orthonormal_vec(vec1);
                Vec3 vec3 = (vec1.cross(vec2)).normalized();
                U(0, dim_b) = vec2[0];
                U(1, dim_b) = vec2[1];
                U(2, dim_b) = vec2[2];
                U(0, dim_c) = vec3[0];
                U(1, dim_c) = vec3[1];
                U(2, dim_c) = vec3[2];
                if (U.det() < 0.0)
                {
                    U(0, dim_b) *= -1.0;
                    U(1, dim_b) *= -1.0;
                    U(2, dim_b) *= -1.0;
                }
                done = 1;
                break; // out of for
            }
        }

        // handle the case where one singular value is small,
        // but the other two are not
        // handle it by computing the cross product of the two eigenvectors
        // for the two large singular values
        if(!done)
        {
            for(int dim = 0; dim < 3; dim++)
            {
                int dim_a = dim;
                int dim_b = (dim + 1) % 3;
                int dim_c = (dim + 2) % 3;

                if(sigma[dim_a] < eps)
                {
                    // columns dimB and dimC are both good,
                    // but column dimA corresponds to a tiny singular value

                    Vec3 vec1(U(0,dim_b), U(1,dim_b), U(2,dim_b)); // column dimB
                    Vec3 vec2(U(0,dim_c), U(1,dim_c), U(2,dim_c)); // column dimC
                    Vec3 vec3 = (vec1.cross(vec2)).normalized();
                    U(0, dim_a) = vec3[0];
                    U(1, dim_a) = vec3[1];
                    U(2, dim_a) = vec3[2];
                    if(U.det() < 0.0)
                    {
                        U(0, dim_a) *= -1.0;
                        U(1, dim_a) *= -1.0;
                        U(2, dim_a) *= -1.0;
                    }
                    done = 1;
                    break; // out of for
                }
            }
        }

        if ((!done) && (modifiedSVD == 1))
        {
            // Handle situation:
            // 3. negative determinant (Tet is inverted in solid mechanics)
            //    - check if det(U) == -1
            //    - If yes, then negate the minimal element of Sigma
            //      and the corresponding column of U

            float det_U = U.det();
            if (det_U < 0.0)
            {
                // negative determinant
                // find the smallest singular value (they are all non-negative)
                int smallest_singular_value_idx = 0;
                for(int dim=1; dim<3; dim++)
                    if (sigma[dim] < sigma[smallest_singular_value_idx])
                        smallest_singular_value_idx = dim;

                // negate the smallest singular value
                sigma[smallest_singular_value_idx] *= -1.0;
                U(0, smallest_singular_value_idx) *= -1.0;
                U(1, smallest_singular_value_idx) *= -1.0;
                U(2, smallest_singular_value_idx) *= -1.0;
            }
        }
    }

    return 0;
}
