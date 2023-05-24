#include "hpc.h"
#include <mpi.h>

void mpi_sum_crosspoints(coupling_data* coupling, double* x, double* cp_buffer) {
    index n_cp = coupling -> n_local_cp;
    index* cp = coupling -> crossPts;
    index* l2g = coupling -> l2g;
    for (index i = 0; i < n_cp; i++) {
        cp_buffer[l2g[cp[i]]] = x[cp[i]];
    }
    MPI_Allreduce(MPI_IN_PLACE, cp_buffer, coupling->n_global_cp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (index i = 0; i < n_cp; i++) {
        x[cp[i]] = cp_buffer[l2g[cp[i]]];
    }
}

double mpi_dotprod(index n, double* x, double* y) {
    double dotprod = 0;
    for (int i = 0; i < n; i++) {
        dotprod += x[i]*y[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotprod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return dotprod;
}