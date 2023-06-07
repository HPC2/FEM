#include "hpc.h"

void mpi_gs(sed* A_sparse, coupling_data* coupling, comm_buffers* buffers, mesh* local_mesh, double* x, double* b) {
    
    index n = A_sparse->n;
    index* fixed = local_mesh->fixed;
    index nfixed = local_mesh->nfixed;
    double tol = 1e-8;
    gem* A = sed_to_dense(A_sparse, true);
    index incRowA = A->incRow;
    index incColA = incColA;
    double* A_data = A->x;

    index* v_nodes = coupling->crossPts;
    index* e_nodes = coupling->e_nodes;
    index* i_nodes = coupling->i_nodes;
    index n_v_nodes = coupling->n_local_cp;
    index n_e_nodes = coupling->n_e_nodes;
    index n_i_nodes = coupling->n_i_nodes;
    
    double omega = 0.2;

    double* resi = calloc(n, sizeof(double));       // get workspace for residual
    double* w    = calloc(n, sizeof(double));       // get temporary workspace
    double* s    = calloc(n, sizeof(double));       // get temporary workspace

    

    double* d = accumulate_inv_diag(coupling, A_sparse, buffers);
    dscal(n, omega, d, 1);
    // resi <-- b
    dcopy(n, b, 1, resi, 1);
    // resi <-- b(aka resi) - A*x
    sysed_spmv(-1, A_sparse, x, 1, 1, resi, 1);
    // set fixed nodes to zero
    for ( size_t i = 0; i < nfixed; ++i){
        resi[fixed[i]] = 0;
    }
    dcopy(n, resi, 1, w, 1); // w <- resi
    mpi_convert_type2_to_type1(coupling, w, buffers);  // convert w
    double sigma = mpi_dotprod(n, w, resi);
    double sigma0 = sigma;


    while (sigma > tol*sigma0) {
        dcopy(n, b, 1, resi, 1);
        
        // v nodes
        indexed_dgemv(n_v_nodes, n_v_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, v_nodes, v_nodes);
        indexed_dgemv(n_v_nodes, n_e_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, e_nodes, v_nodes);
        indexed_dgemv(n_v_nodes, n_i_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, i_nodes, v_nodes);

        indexed_dcopy(n_v_nodes, resi, 1, w, 1, v_nodes);
        
        mpi_sum_crosspoints(coupling, w, buffers->cp_buffer);
        
        indexed_dcopy(n_v_nodes, w, 1, s, 1, v_nodes);
        indexed_dmult(n_v_nodes, d, 1, s, 1, v_nodes);
        
        for ( index i = 0; i < nfixed; ++i){
            s[fixed[i]] = 0;
        }
        return;
        indexed_daxpy(n_v_nodes, 1, s, 1, x, 1, v_nodes);

        

        // e nodes
        indexed_dgemv(n_e_nodes, n_v_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, v_nodes, e_nodes);
        indexed_dgemv(n_e_nodes, n_e_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, e_nodes, e_nodes);
        indexed_dgemv(n_e_nodes, n_i_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, i_nodes, e_nodes);
        indexed_dcopy(n_e_nodes, resi, 1, w, 1, e_nodes);
        mpi_sum_crosspoints(coupling, w, buffers->cp_buffer);
        indexed_dcopy(n_e_nodes, w, 1, s, 1, e_nodes);
        indexed_dmult(n_e_nodes, d, 1, s, 1, e_nodes);
        for ( size_t i = 0; i < nfixed; ++i) {
            s[fixed[i]] = 0;
        }
        indexed_daxpy(n_e_nodes, 1, s, 1, x, 1, e_nodes);


        // i nodes 
        indexed_dgemv(n_i_nodes, n_v_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, v_nodes, i_nodes);
        indexed_dgemv(n_i_nodes, n_e_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, resi, 1, e_nodes, i_nodes);


        // gauss-seidel step
        
        

        for (index i = 0; i < n_i_nodes; i++) {
            
            double dx = resi[i_nodes[i]];
            
            for (index j = 0; j < n_i_nodes; j++) {
                if (i != j)
                    dx -= A_data[i_nodes[i]*incRowA + i_nodes[j]*incColA] * x[i_nodes[j]];
            }
            x[i_nodes[i]] += d[i_nodes[i]] * dx;
        }

        

        // end gauss-seidel
        indexed_dcopy(n_i_nodes, resi, 1, w, 1, i_nodes);
        indexed_dgemv(n_i_nodes, n_i_nodes,
            -1, A_data, incRowA, incColA,
            x, 1,
            1, w, 1, i_nodes, i_nodes);
        sigma = mpi_dotprod(n, w, resi);
    }
}