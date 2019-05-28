#include "eigen_3x3.hpp"

#include <cmath>

const int n = 3;

// -----------------------------------------------------------------------------

static float hypot2(float x, float y) {
    return sqrtf(x*x+y*y);
}

// -----------------------------------------------------------------------------

// Symmetric Householder reduction to tridiagonal form.
static void tred2(float V[3][3], float d[3], float e[3]) {

    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int j = 0; j < n; j++) {
        d[j] = V[n-1][j];
    }

    // Householder reduction to tridiagonal form.

    for (int i = n-1; i > 0; i--) {

        // Scale to avoid under/overflow.

        float scale = 0.0;
        float h = 0.0;
        for (int k = 0; k < i; k++) {
            scale = scale + fabs(d[k]);
        }
        if (scale == 0.0) {
            e[i] = d[i-1];
            for (int j = 0; j < i; j++) {
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        } else {

            // Generate Householder vector.

            for (int k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            float f = d[i-1];
            float g = sqrtf(h);
            if (f > 0) {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (int j = 0; j < i; j++) {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for (int j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (int k = j+1; k <= i-1; k++) {
                    g += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            float hh = f / (h + h);
            for (int j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }
            for (int j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i-1; k++) {
                    V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n-1; i++) {
        V[n-1][i] = V[i][i];
        V[i][i] = 1.0;
        float h = d[i+1];
        if (h != 0.0) {
            for (int k = 0; k <= i; k++) {
                d[k] = V[k][i+1] / h;
            }
            for (int j = 0; j <= i; j++) {
                float g = 0.0;
                for (int k = 0; k <= i; k++) {
                    g += V[k][i+1] * V[k][j];
                }
                for (int k = 0; k <= i; k++) {
                    V[k][j] -= g * d[k];
                }
            }
        }
        for (int k = 0; k <= i; k++) {
            V[k][i+1] = 0.0;
        }
    }
    for (int j = 0; j < n; j++) {
        d[j] = V[n-1][j];
        V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
}

// -----------------------------------------------------------------------------

// Symmetric tridiagonal QL algorithm.
static void tql2(float V[3][3], float d[3], float e[3]) {

    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int i = 1; i < n; i++) {
        e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    float f = 0.0f;
    float tst1 = 0.0f;
    float eps = (float)pow(2.0,-52.0);
    for (int l = 0; l < n; l++) {

        // Find small subdiagonal element

        tst1 = fmax(tst1, fabs(d[l]) + fabs(e[l]));
        int m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps*tst1) {
                break;
            }
            m++;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.

        if (m > l) {
            int iter = 0;
            do {
                iter = iter + 1;  // (Could check iteration count here.)

                // Compute implicit shift

                float g = d[l];
                float p = (d[l+1] - g) / (2.0f * e[l]);
                float r = hypot2(p, 1.0f);
                if (p < 0) {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                float dl1 = d[l+1];
                float h = g - d[l];
                for (int i = l+2; i < n; i++) {
                    d[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d[m];
                float c = 1.0;
                float c2 = c;
                float c3 = c;
                float el1 = e[l+1];
                float s = 0.0;
                float s2 = 0.0;
                for (int i = m-1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p,e[i]);
                    e[i+1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i+1] = h + s * (c * g + s * d[i]);

                    // Accumulate transformation.

                    for (int k = 0; k < n; k++) {
                        h = V[k][i+1];
                        V[k][i+1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence.

            } while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for (int i = 0; i < n-1; i++) {
        int k = i;
        float p = d[i];
        for (int j = i+1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}

// -----------------------------------------------------------------------------

void eigen_decomposition(const float A[3][3], float V[3][3], float d[3]) {
    float e[n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V[i][j] = A[i][j];
        }
    }
    tred2(V, d, e);
    tql2(V, d, e);
}
