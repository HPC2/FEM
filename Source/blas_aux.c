#include "hpc.h"
#include <math.h>           // for nan(), fabs()
#include <stdlib.h>         // for malloc(), free()

void
printfDGeMatrix(const char * fmt, size_t m, size_t n,
                const double *A,
                ptrdiff_t incRowA, ptrdiff_t incColA)
{
    for (size_t i=0; i<m; ++i) {
        for (size_t j=0; j<n; ++j) {
            printf(fmt, A[i*incRowA + j*incColA]);
        }
        printf("\n");
    }
    printf("\n");
}

void
printDGeMatrix(size_t m, size_t n,
               const double *A,
               ptrdiff_t incRowA, ptrdiff_t incColA)
{
    printfDGeMatrix("%9.4lf ", m, n, A, incRowA, incColA);
}

void
printIGeMatrix(size_t m, size_t n,
               const size_t *A,
               ptrdiff_t incRowA, ptrdiff_t incColA)
{
    for (size_t i=0; i<m; ++i) {
        for (size_t j=0; j<n; ++j) {
            printf("%6zu ", A[i*incRowA + j*incColA]);
        }
        printf("\n");
    }
    printf("\n");
}


