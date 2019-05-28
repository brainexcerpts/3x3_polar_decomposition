#include "polar_decomposition_3x3.hpp"

#include "svd_3x3.hpp"
#include "vec3.hpp"

#include <cstdio>

// =============================================================================
namespace details {
// =============================================================================

/// a, b, c are 3d-vectors
/// compute cross product c = a x b
static inline
void cross_product(const float* a, const float* b, float* c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

// -----------------------------------------------------------------------------

/// one-norm of a 3 x 3 matrix (row major)
static inline
float one_norm(const float* A)
{
    float norm = 0.0f;
    for (int i = 0; i < 3; i++)
    {
        float col_abs_sum = fabs(A[i + 0]) + fabs(A[i + 3]) + fabs(A[i + 6]);
        if (col_abs_sum > norm)
            norm = col_abs_sum;
    }
    return norm;
}

// -----------------------------------------------------------------------------

/// infinity-norm of a 3 x 3 matrix (row major)
static inline
float inf_norm(const float* A)
{
    float norm = 0.0;
    for (int i = 0; i < 3; i++)
    {
        float row_sum = fabs(A[3 * i + 0]) + fabs(A[3 * i + 1]) + fabs(A[3 * i + 2]);
        if (row_sum > norm)
            norm = row_sum;
    }
    return norm;
}

} // Namespace details =========================================================


// Input: M (3x3 mtx)
// Note All matrices are row-major in this implementation
template<bool force_rotation>
float Polar_decomposition<force_rotation>::compute(const Mat3& M, float tolerance)
{
    using namespace details;
    Mat3 R;
    Mat3 S;
    //------

    float Mk[9];
    float Ek[9];
    float det, M_one_norm, M_inf_norm, E_one_norm;
    bool use_svd = false;

    // Mk = M^T
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            Mk[3 * i + j] = M[3 * j + i];

    M_one_norm = one_norm(Mk);
    M_inf_norm = inf_norm(Mk);

    do
    {
        float M_adj_Tk[9];

        // row 2 x row 3
        cross_product(&(Mk[3]), &(Mk[6]), &(M_adj_Tk[0]));
        // row 3 x row 1
        cross_product(&(Mk[6]), &(Mk[0]), &(M_adj_Tk[3]));
        // row 1 x row 2
        cross_product(&(Mk[0]), &(Mk[3]), &(M_adj_Tk[6]));

        det = Mk[0] * M_adj_Tk[0] + Mk[1] * M_adj_Tk[1] + Mk[2] * M_adj_Tk[2];

        if( det <= 1e-6 )
        {
            use_svd = true;
            break;
        }

        if( det == 0.0 )
        {
            printf("Warning (polarDecomposition) : zero determinant encountered.\n");
            break;
        }

        float MadjT_one = one_norm(M_adj_Tk);
        float MadjT_inf = inf_norm(M_adj_Tk);

        float gamma = sqrt(sqrt((MadjT_one * MadjT_inf) / (M_one_norm * M_inf_norm * det * det)));
        float g1 = gamma * 0.5f;
        float g2 = 0.5f / (gamma * det);

        for(int i = 0; i < 9; i++)
        {
            Ek[i] = Mk[i];
            Mk[i] = g1 * Mk[i] + g2 * M_adj_Tk[i];
            Ek[i] -= Mk[i];
        }

        E_one_norm = one_norm(Ek);
        M_one_norm = one_norm(Mk);
        M_inf_norm = inf_norm(Mk);

    } while ( E_one_norm > M_one_norm * tolerance );

    if(use_svd && force_rotation)
    {
        // use the SVD algorithm to compute R
        Mat3 Mm = M;
        Mat3 Um, Vm;
        Vec3 lambda;
        int modified_SVD = 1;
        SVD_3x3(Mm, Um, lambda, Vm, tolerance, modified_SVD);
        R = Um * Vm.transpose();
    }
    else
    {
        // R = Mk^T
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                R[3*i + j] = Mk[3*j + i];
    }

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            S[3*i + j] = 0;
            for(int k = 0; k < 3; k++)
                S[3*i + j] += Mk[3*i + k] * M[3*k + j];
        }

    // S must be symmetric; enforce the symmetry
    S[1] = S[3] = 0.5f * (S[1] + S[3]);
    S[2] = S[6] = 0.5f * (S[2] + S[6]);
    S[5] = S[7] = 0.5f * (S[5] + S[7]);

    // --------
    _R = R;
    _S = S;
    return det;
}
