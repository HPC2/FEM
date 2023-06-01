#include "hpc.h"
#include <mpi.h>

void mpi_cg(sed* A, coupling_data* coupling, comm_buffers* buffers, mesh* local_mesh, double* x, double* b) {
      int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);

      index n = A->n;
      index* fixed = local_mesh->fixed;
      index nfixed = local_mesh->nfixed;
      double tol = 1e-6;

      double* resi = calloc (n, sizeof(double));       // get workspace for residual
      double* w    = calloc (n, sizeof(double));       // get temporary workspace
      double* s    = calloc (n, sizeof(double));       // get temporary workspace
      double* v    = calloc (n, sizeof(double));       // get temporary workspace

      // cg solver

      double* inv_diag = accumulate_inv_diag(coupling, A, buffers);  
      dcopy(n, b, 1, resi, 1); // resi <- b
      sysed_spmv(-1, A, x, 1, 1, resi, 1);  // resi <- b - A*x
      for ( size_t i = 0; i < nfixed; ++i){
            resi[fixed[i]] = 0;
      }
      dcopy(n, resi, 1, w, 1); // w <- resi
      mpi_convert_type2_to_type1(coupling, w, buffers);  // convert w
      // dmult(n, inv_diag, 1, w, 1); // conditioning: w <- D^-1 * w
      dcopy(n, w, 1, s, 1); // copy w to s
      double sigma = mpi_dotprod(n, w, resi);
      double sigma0 = sigma;
      double sigma_old = sigma;
      double alpha, beta;
      index iter = 0;
      while (sigma > tol*sigma0) {
        sysed_spmv(1, A, s, 1, 0, v, 1);  // v <- A*s
        for ( size_t i = 0; i < nfixed; ++i){
              v[fixed[i]] = 0;
        }
        alpha = sigma/mpi_dotprod(n, s, v); // alpha -> sigma/(s, v)
        daxpy(n, alpha, s, 1, x, 1); // x <- alpha*s + x
        daxpy(n, -alpha, v, 1, resi, 1); // r <- -alpha*v + r
        // recalc w
        dcopy(n, resi, 1, w, 1); // w <- resi
        mpi_convert_type2_to_type1(coupling, w, buffers);  // convert w
        // dmult(n, inv_diag, 1, w, 1); // conditioning: w <- D^-1 * w
        // recalc error 
        sigma = mpi_dotprod(n, w, resi);
        //if (rank == 0) printf("sigma: %4.4lf\n", sigma);
        beta = sigma/sigma_old;
        //if (rank == 0) printf("beta: %4.4lf\n", beta);
        sigma_old = sigma;
        dscal(n, beta, s, 1);
        daxpy(n, 1, w, 1, s, 1);
      }

      free(inv_diag);
      free(resi);
      free(w);
      free(s);
      free(v);
}