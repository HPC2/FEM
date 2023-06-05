#include "hpc.h"

index seq_gs(sed* A, mesh* Mesh, double* x, double* b) {

    index n = A->n;
    double* w    = calloc (n, sizeof(double));
    double* resi = calloc (n, sizeof(double));
    index nfixed = Mesh->nfixed ; 
    index* fixed       = Mesh->fixed ;

    double tol = 1e-8;

    index n_iter = 0;

    // Calculate square norm of residual
     // resi <-- b
     dcopy(n, b, 1, resi, 1);
     // resi <-- b - A*x
     sysed_spmv(-1, A, x, 1, 1, resi, 1);
     // set fixed nodes to zero
     for ( size_t i = 0; i < nfixed; ++i){
         resi[fixed[i]] = 0;
     }
     // resi_norm <-- || resi ||_2^2
     double sigma = ddot(n, resi, 1, resi, 1);
     double sigma0 = sigma;

    while (sigma > tol*sigma0) {
        n_iter++;
        // Sym Gauss-Seidel iterations
        sed_gs_constr(A, b, x, w, fixed, nfixed, 1); 
        sed_gs_constr(A, b, x, w, fixed, nfixed, 0); 
        // Calculate square norm of residual
        // resi <-- b
        dcopy(n, b, 1, resi, 1);
        // resi <-- b - A*x
        sysed_spmv(-1, A, x, 1, 1, resi, 1);
        // set fixed nodes to zero
        for ( size_t i = 0; i < nfixed; ++i){
            resi[fixed[i]] = 0;
        }
        // resi_norm <-- || resi ||_2^2
        sigma = ddot(n, resi, 1, resi, 1);
    }

    free(w);
    free(resi);
    return n_iter;
}