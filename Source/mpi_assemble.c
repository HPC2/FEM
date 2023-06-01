#include "hpc.h"
#include <mpi.h>

double* mpi_assemble_t2_vec(coupling_data* coupling, double* local_x, index n) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    index* l2g = coupling->l2g;
    double* global_x = calloc(coupling->n_global_nodes, sizeof(double));
    for (index i = 0; i < n; i++) {
        global_x[l2g[i]] = local_x[i];
    }
    if (rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, global_x, coupling->n_global_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        return global_x;
    } else {
        MPI_Reduce(global_x, global_x, coupling->n_global_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        free(global_x);
        return NULL;
    }
}

double* mpi_assemble_t1_vec(coupling_data* coupling, double* local_x, index n) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    index* l2g = coupling->l2g;
    double* global_x = calloc(coupling->n_global_nodes, sizeof(double));
    
    for (index i = 0; i < n; i++) {
        global_x[l2g[i]] = local_x[i];
    }
    if (rank == 0) {
        double* max_x = calloc(coupling->n_global_nodes, sizeof(double));
        double* min_x = calloc(coupling->n_global_nodes, sizeof(double));
        MPI_Reduce(global_x, max_x, coupling->n_global_nodes, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(global_x, min_x, coupling->n_global_nodes, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        for (int i = 0; i < coupling->n_global_nodes; i++) {
            if (fabs(max_x[i] - min_x[i]) < 1e-6) {
                global_x[i] = max_x[i];
            } else {
                global_x[i] = max_x[i] + min_x[i];
            }
        }
        free(max_x);
        free(min_x);
        return global_x;
    } else {
        MPI_Reduce(global_x, global_x, coupling->n_global_nodes, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(global_x, global_x, coupling->n_global_nodes, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        free(global_x);
        return NULL;
    }
}

double* mpi_assemble_A (sed* A_loc, coupling_data* coupling){
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    index* l2g = coupling->l2g;
    index n_loc = A_loc->n;
    index n_glob = coupling->n_global_nodes;
    double* A_global = calloc(n_glob*n_glob, sizeof(double)); // rowmajor
    gem* A_loc_full = sed_to_dense (A_loc, true);
    double* A_loc_data = A_loc_full->x;
    index incRow = A_loc_full->incRow;
    index incCol = A_loc_full->incCol;
   
    //index* l2g = coupling->l2g;
    for (index i = 0; i<n_loc; i++){
        for (index j=0; j<n_loc; j++)
            A_global[l2g[i]*n_glob + l2g[j]] = A_loc_data[i*incRow + j*incCol];
    }
    
    if (rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, A_global, n_glob*n_glob, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        return A_global;
    } else {
        MPI_Reduce(A_global, A_global, n_glob*n_glob , MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        gem_free(A_loc_full);
        free(A_global);
        return NULL;
    }
}
