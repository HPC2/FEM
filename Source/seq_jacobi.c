#include "hpc.h"

void seq_jacobi(sed* A, mesh* Mesh, double* x, double* b) {

    index n = A->n;
    double* w    = calloc (n, sizeof(double));
    double* resi = calloc (n, sizeof(double));
    double* d = malloc((n)*sizeof(double));
    index nfixed = Mesh->nfixed ; 
    index* fixed = Mesh->fixed ;

    // Jacobi sequetiell
    //--------------------------------------------------------//
    // choose tol
    double tol = 1e-8;
    // choose omega 
    double omega = 0.2;
    //--------------------------------------------------------//
    // resi <-- b
    dcopy(n, b, 1, resi, 1);
    // resi <-- b(aka resi) - A*x
    sysed_spmv(-1, A, x, 1, 1, resi, 1);
    // set fixed nodes to zero
    for ( size_t i = 0; i < nfixed; ++i){
        resi[fixed[i]] = 0;
    }
    // d = diag(D)^-1
    double* A_data = A->x;
    for(int i=0; i<n; i++){
        d[i] = 1/A_data[i]; 
    }
    // w <- resi
    dcopy(n, resi,1, w,1);
    // w <- d.*resi(aka w)
    dmult(n, d,1, w,1);
    // sigma <- w'*resi
    double sigma0 = ddot(n, w,1, resi,1);
    double sigma = sigma0;
    int k = 0;
    //for(k = 0; k<10000; ) {
    while(sigma > tol*sigma0){
        k = k+1;
        // x <- x + omega*w
        daxpy(n, omega,w,1, x,1);
        // resi <-- b
        dcopy(n, b, 1, resi, 1);
        // resi <-- b(aka resi) - A*x
        sysed_spmv(-1,A, x,1, 1, resi,1);
        // set fixed nodes to zero
        for ( size_t i = 0; i < nfixed; ++i){
            resi[fixed[i]] = 0;
        }
        // w <- resi
        dcopy(n, resi,1, w,1);
        // w <- d.*resi(aka w)
        dmult(n, d,1, w,1);
        // sigma <- w'*resi
        sigma = ddot(n, w,1, resi,1);
    }
    // end Jacobi sequentiell

    free(w);
    free(resi);
    free(d);
}