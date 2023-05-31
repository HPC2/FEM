#include "hpc.h"

void cg_seq(const sed * A, double* f, double *u, const double tol, double* u_D, const mesh *M) {    
    index n =  A -> n;
    double r[n];
    double v[n];
    double alpha;
    double beta;
    
    double w[n];
    double s[n];
/*
    set_dir(M,u);
    sysed_spmv(-1,A,u_D,1,1,f,1);
    dcopy(n,r,1,f,1);
    sysed_spmv(-1,A,u,1,1,r,1);
    set_dir(M,r);
*/
    dcopy(n,f,1,r,1);
    sysed_spmv(-1,A,u,1,1,r,1);

    double* C = malloc((n)*sizeof(double));
    dcopy(n, A->x, 1, C, 1);
    for(int i=0; i<n; i++){
        C[i] = 1/C[i];
    }
    dcopy(n, r,1, w,1); // w <- resi
    dmult(n, C,1, w,1); // w <- d.*w = d.*resi
    

    dcopy(n, w,1, s,1);

    double sigma_0 = ddot(n,w,1,r,1);
    double sigma = sigma_0;
    double sigma_old = sigma_0;
    double dd;
    
    while ( sqrt(sigma/sigma_0) > tol ) {

        printf("segma/sigma_0 = %lf\n",sigma/sigma_0);
        sysed_spmv(1,A,s,1,1,v,1);
        
        dd = ddot(n,s,1,v,1);
        set_dir(M,v);
        alpha = sigma / dd;
        
        daxpy(n,alpha,s,1,u,1);
        daxpy(n,alpha*(-1),v,1,r,1);
        dcopy(n, r,1, w,1); // w <- resi
        dmult(n, C,1, w,1);

        sigma = ddot(n,w,1,r,1);

        beta = sigma/sigma_old;
        sigma_old = sigma;

        dscal(n,beta,s,1);
        daxpy(n,1,w,1,s,1);
        
    }  

    /*set_dir_u_D(M,u,u_D);*/
} 













