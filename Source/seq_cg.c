#include "hpc.h"

index seq_cg(sed* A, mesh* Mesh, double* x, double* b) {

    index n = A->n;
    double* w    = calloc (n, sizeof(double));
    double* resi = calloc (n, sizeof(double));
    double* s    = calloc (n, sizeof(double));       // get temporary workspace
    double* v    = calloc (n, sizeof(double));       // get temporary workspace
    index nfixed = Mesh->nfixed ; 
    index* fixed = Mesh->fixed ;

    index n_iter = 0;

    double tol = 1e-8;

    dcopy(n, b, 1, resi, 1);
    sysed_spmv(-1, A, x, 1, 1, resi, 1);
    for ( size_t i = 0; i < nfixed; ++i){
            resi[fixed[i]] = 0;
    }
    dcopy(n, resi, 1, w, 1);
    dcopy(n, w, 1, s, 1);
    double sigma = ddot(n, w, 1, resi, 1);
    double sigma0 = sigma;
    double sigma_old = sigma;
    double alpha;
    double beta;
    while (sigma > tol*sigma0) {
        n_iter++;
        sysed_spmv(1, A, s, 1, 0, v, 1);
        for ( size_t i = 0; i < nfixed; ++i){
            v[fixed[i]] = 0;
        }
        alpha = sigma / ddot(n, s, 1, v, 1);
        daxpy(n, alpha, s, 1, x, 1);
        daxpy(n, -alpha, v, 1, resi, 1);
        dcopy(n, resi, 1, w, 1);
        sigma = ddot(n, w, 1, resi, 1);
        // printf("sigma: %1.4lf\n", sigma);
        beta = sigma/sigma_old;
        sigma_old = sigma;
        dscal(n, beta, s, 1);
        daxpy(n, 1, w, 1, s, 1);
    }

    free(w);
    free(resi);
    free(v);
    free(s);
    return n_iter;
}