#include "hpc.h"
#include <math.h>           // for nan(), fabs()
#include <stdlib.h>         // for malloc(), free()

void
dcopy(index n,
      const double *x, index incX,
      double       *y, index incY)
{
    for (index i=0; i<n; ++i) {
        y[i*incY] = x[i*incX];
    }
}

void
daxpy(index n, double alpha,
      const double *x, index incX,
      double       *y, index incY)
{
    if (alpha == 0){
        return;
    }
    for (index i=0; i<n; ++i) {
        y[i*incY] += alpha*x[i*incX];
    }
}

double
ddot(index n,
     const double *x, index incX,
     const double *y, index incY)
{
    double result = 0;

    for (index i=0; i<n; ++i) {
        result += x[i*incX]*y[i*incY];
    }

    return result;
}

index
idamax(index n, const double *x, index incX)
{
    double max = 0;
    index I    = 0;

    for (index i=0; i<n; ++i) {
        if (fabs(x[i*incX])>max) {
            I   = i;
            max = fabs(x[i*incX]);
        }
    }
    return I;
}

void
dswap(index n, 
      double *x, index incX, 
      double *y, index incY)
{
    for (index i=0; i<n; ++i) {
        double tmp = x[i*incX];
        x[i*incX] = y[i*incY];
        y[i*incY] = tmp;
    }
}

// x <- alpha * x
void
dscal(index  n,
      double alpha,
      double *x, index incX)
{
    // No-op if alpha 1
    if (alpha==1) {
        return;
    }

    // Scaling
    if (alpha!=0) {
        for (index i=0; i<n; ++i) {
            x[i*incX] *= alpha;
        }
    } else {
        for (index i=0; i<n; ++i) {
            x[i*incX] = 0;
        }
    }
}


