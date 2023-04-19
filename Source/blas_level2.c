#include "hpc.h"

// symmetric sed general matrix vector product
// y <- alpha * A * x + beta * y 
void
sysed_spmv(double alpha,
           const sed *A,
           const double *x, index incX,
           double beta,
           double *y, index incY)
{
    index p, j, n, *Ai ;
    double *Ax;

    // For convenience
    n  = A->n ; 
    Ai = A->i ; 
    Ax = A->x ;

    // Scal y
    dscal(n, beta, y, incY);

    if (alpha==0 || n==0) {
        return;
    }

    for (j = 0 ; j < n ; j++)
    {
        // Diagonal
        y[j*incY] += alpha*Ax[j] * x[j*incX] ;

        // Sub/Super-Diagonal
        for (p = Ai[j] ; p < Ai[j+1] ; p++)
        {
            // lower part axpy-based
            y[Ai[p]*incY] += alpha * Ax[p] * x[j*incX] ;

            // upper part ddot-based
            y[j*incY]     += alpha * Ax[p] * x[Ai[p]*incX] ;
        }
    }
}


#ifndef AXPYF
#define AXPYF 4
#endif

// y <- alpha*A*x + y where A is a m x AXPYF matrix
void
daxpyf(size_t m, double alpha,
       const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
       const double *x, ptrdiff_t incX,
       double *y, ptrdiff_t incY)
{
    for (size_t i=0; i<m; ++i) {
        for (size_t l=0; l<AXPYF; ++l) {
            y[i*incY] += alpha * A[i*incRowA+l*incColA] * x[l*incX];
        }
    }
}

void
dgemv_axpyf(size_t m, size_t n,
            double alpha,
            const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
            const double *x, ptrdiff_t incX,
            double *y, ptrdiff_t incY)
{
    size_t nb = n / AXPYF;
    for (size_t j=0; j<nb; ++j) {
        daxpyf(m, alpha,
               &A[j*AXPYF*incColA], incRowA, incColA,
               &x[j*AXPYF*incX], incX,
               y, incY);
    }
    for (size_t j=nb*AXPYF; j<n; ++j) {
        daxpy(m, alpha*x[j*incX], &A[j*incColA], incRowA, y, incY);
    }
}

#ifndef DOTF
#define DOTF 4
#endif

void
dgemv_dotf(size_t m, size_t n,
           double alpha,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *x, ptrdiff_t incX,
           double *y, ptrdiff_t incY)
{
    size_t mb = m / DOTF;

    for (size_t i=0; i<mb; ++i) {
        for (size_t j=0; j<n; ++j) {
            for (size_t l=0; l<DOTF; ++l) {
                y[(DOTF*i+l)*incY]
                    += alpha*A[(DOTF*i+l)*incRowA+j*incColA]*x[j*incX];
            }
        }
    }

    for (size_t i=mb*DOTF; i<m; ++i) {
        y[i*incY] += alpha*ddot(n, &A[i*incRowA], incColA, x, incX);
    }
}

void
dgemv(size_t m, size_t n,
      double alpha,
      const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
      const double *x, ptrdiff_t incX,
      double beta,
      double *y, ptrdiff_t incY)
{
    dscal(m, beta, y, incY);

    if (alpha==0 || m==0 || n==0) {
        return;
    }

    if (incRowA<incColA) {
        dgemv_axpyf(m, n, alpha, A, incRowA, incColA, x, incX, y, incY);
    } else {
        dgemv_dotf(m, n, alpha, A, incRowA, incColA, x, incX, y, incY);
    }
}


void
dtrsv(size_t n, bool lower, bool unit,
      const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
      double *x, ptrdiff_t incX)
{
    if (lower) {
        if (incRowA<incColA) {
            // A is col major

            size_t nb = n / AXPYF;
            for (size_t j=0; j<nb; ++j) {
                for (size_t l=0; l<AXPYF; ++l) {
                    size_t J = j*AXPYF + l;

                    if (!unit) {
                        x[J*incX] /= A[J*(incRowA+incColA)];
                    }
                    daxpy(AXPYF-1-l, -x[J*incX],
                          &A[(J+1)*incRowA + J*incColA], incRowA,
                          &x[(J+1)*incX], incX);
                }
                daxpyf(n-(j+1)*AXPYF, -1,
                       &A[((j+1)*incRowA + j*incColA)*AXPYF], incRowA, incColA,
                       &x[j*incX*AXPYF], incX,
                       &x[(j+1)*incX*AXPYF], incX);
            }
            for (size_t j=nb*AXPYF; j<n; ++j) {
                if (!unit) {
                    x[j*incX] /= A[j*(incRowA+incColA)];
                }
                daxpy(n-1-j, -x[j*incX],
                      &A[(j+1)*incRowA+j*incColA], incRowA,
                      &x[(j+1)*incX], incX);
            }
        } else {
            // A is row major
            for (size_t j=0; j<n; ++j) {
                x[j*incX] -= ddot(j, &A[j*incRowA], incColA, x, incX);
                if (!unit) {
                    x[j*incX] /= A[j*(incRowA+incColA)];
                }
            }
        }
    } else {
        if (incRowA<incColA) {
            // A is col major
            for (size_t j=n; j-- >0; ) {
                if (!unit) {
                    x[j*incX] /= A[j*(incRowA+incColA)];
                }
                daxpy(j, -x[j*incX], &A[j*incColA], incRowA,
                      x, incX);
            }
        } else {
            // A is row major
            for (size_t j=n; j-- >0; ) {
                x[j*incX] -= ddot(n-1-j,
                                  &A[j*incRowA+(j+1)*incColA], incColA,
                                  &x[(j+1)*incX], incX);
                if (!unit) {
                    x[j*incX] /= A[j*(incRowA+incColA)];
                }
            }
        }
    }
}
