#include "hpc.h"

void cg_seq(const sed * A, const double* f, double *u, double tol
            const double *u_D, const mesh *M) {    
    index n =  A -> n;
    double r[n];
    double v[n];
    double alpha;
    double beta;
    set_dir(M,u);

    sysed_spmv(-1,A,u_D,f,1);
    dcopy(n,r,1,f,1);
    sysed_spmv(-1,A,u,r,1);
    set_dir(M,r);
    
    double sigma_0 = ddot(n,r,1,r,1);
    double sigma = sigma_0;
    
    while (sqrt(sigma) < tol ) {
        sysed_spmv(-1,A,r,v,1);
        alpha = sigma /ddot(n,r,1,v,1);
        daxpy(n,alpha,r,1,u,1);
        daxpy(n,alpha*(-1),v,1,u,1);
        sigma = ddot(n,r,1,r,1);
        set_dir(M,r);
        beta = sigma/sigma_0;
        sigma = sigma_0;
        daxpy(n,alpha,r,1,r,1);
    }  
    set_dir_u_D(m,u,u_D)
} 















}