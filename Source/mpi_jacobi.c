#include "hpc.h"

void mpi_jacobi(sed* A, coupling_data* coupling, comm_buffers* buffers, mesh* local_mesh, double* x, double* b) {
    
    index n = A->n;
    index* fixed = local_mesh->fixed;
    index nfixed = local_mesh->nfixed;
    double tol = 1e-6;

    double* resi = calloc (n, sizeof(double));       // get workspace for residual
    double* w    = calloc (n, sizeof(double));       // get temporary workspace

    // jacobi solver
    double* inv_diag = accumulate_inv_diag(coupling, A, buffers);
    dcopy(n, b, 1, resi, 1); // resi <- b
    sysed_spmv(-1, A, x, 1, 1, resi, 1);  // resi <- resi - A*x
    dcopy(n, resi, 1, w, 1); // w <- resi
    mpi_convert_type2_to_type1(coupling, w, buffers);  // convert w
    dmult(n, inv_diag, 1, w, 1); // w <- D^-1 * w
    double sigma = mpi_dotprod(n, w, resi);
    double sigma0 = sigma;
    
    double omega = 1; // full jacobi steps for now
    while (sigma > tol*sigma0) {
      // update x
      daxpy(n, omega, w, 1, x, 1);
      // recalc residuum
      dcopy(n, b, 1, resi, 1); // resi <- b
      sysed_spmv(-1, A, x, 1, 1, resi, 1);  // resi <- resi - A*x
      for ( size_t i = 0; i < nfixed; ++i){
            resi[fixed[i]] = 0;
      }
      // recalc w
      dcopy(n, resi, 1, w, 1); // w <- resi
      mpi_convert_type2_to_type1(coupling, w, buffers);  // convert w
      dmult(n, inv_diag, 1, w, 1); // w <- D^-1 * w
      // recalc error 
      sigma = mpi_dotprod(n, w, resi);
    }

    free(inv_diag);
    free(resi);
    free(w);
}