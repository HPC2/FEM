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