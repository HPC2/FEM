#include "hpc.h"

// Gauss-Seidel Iteration with constrains x(fixed) = b(fixed) 
index sed_gs_constr (const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed, bool forward)
{
    // Declare variables
    index p, j, n, t, ft, *Ap, *Ai ;
    double *Ax;

    // check inputs
    if ( !A || !x || !w || !b ) return (0) ;  /* check inputs */

    // For convenience
    n  = A->n ; 
    Ai = A->i ; 
    Ax = A->x ;
    
    // w <-- b
    dcopy(n, b, 1, w, 1);

    // Forward Gauss Iteration
    if (forward){
        // w <-- b - U*x_k
        for (j = 0 ; j < n; j++){
            for (p = Ai[j] ; p < Ai[j+1] ; p++){
                w[j] -= Ax[p] * x [Ai[p]];
            }
        }

        // Initialize fixed bdry variables (from first to last)
        ft = fixed[0]; 
        t  = 1; 

        // If not dirichlet: x_k+1 <-- (L+D)^-1 * (b - U*x_k)
        // Basically a triangular solve
        for (j = 0 ; j < n; j++)
        {
            if ( j != ft){
                x[j] = w[j] / Ax[j];
            } else {
                if (t < nFixed) ft = fixed[t++];
            }
            for (p = Ai[j] ; p < Ai[j+1] ; p++){
                w[Ai[p]] -= Ax[p] * x [j];
            }
        }
    // Backward Gauss Iteration
    } else {
        // w <-- b - L*x_k
        for (j = 0 ; j < n; j++){
            for (p = Ai[j] ; p < Ai[j+1] ; p++){
                w[Ai[p]] -= Ax[p] * x [j];
            }
        }
        // Initialize fixed bdry variables (now from last to first)
        t  = nFixed-1; 
        ft = fixed[nFixed-1]; 
        // If not dirichlet: x_k+1 <-- (D+U)^-1 * (b - L*x_k)
        // Basically a triangular solve
        for (j = n-1 ; j >=0; j--){
            if ( j != ft){
                for (p = Ai[j] ; p < Ai[j+1] ; p++){
                    w[j] -= Ax[p] * x [Ai[p]];
                }
                x[j] = w[j] / Ax[j];
            } else {
                if (t > 0 ){
                    ft = fixed[--t];
                }
            } 
        }
    }    

    return (1) ;
}
