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

void mpi_sum_interfaces(coupling_data* coupling, double* x, double* if_buffer_send, double* if_buffer_recv) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for(int i=0; i<4; i++) { // iterate through all 4 colors
        index left  = coupling->coupl[5*i+2];
        index right = coupling->coupl[5*i+3];
        index color = coupling->coupl[5*i+4];
        index* if_nodes = coupling->icoupl[i]; // if non-exitent, this is NULL
        index n_nodes_interface = coupling->dcoupl[i];
        if(color == -1) {
            left  = MPI_PROC_NULL;
            right = MPI_PROC_NULL;
        }
        
        if(color > -1) {
            for (index j = 0; j < n_nodes_interface; j++) {
                if_buffer_send[j] = x[if_nodes[j]];
            }
        }

        if(rank == left) { // I am left and communicate with right
            MPI_Send(if_buffer_send, n_nodes_interface, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);
            MPI_Recv(if_buffer_recv, n_nodes_interface, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }else{             // I am right and communicate with left
            MPI_Recv(if_buffer_recv, n_nodes_interface, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(if_buffer_send, n_nodes_interface, MPI_DOUBLE, left, 0, MPI_COMM_WORLD);
        }

        if(color > -1) {
            for (index j = 0; j < n_nodes_interface; j++) {
                x[if_nodes[j]] += if_buffer_recv[j]; // accumulating
            }
        }
    }
}

void mpi_convert_type2_to_type1(coupling_data* coupling, double* x, double* cp_buffer, double* if_buffer_send, double* if_buffer_recv) {
    mpi_sum_crosspoints(coupling, x, cp_buffer); 
    mpi_sum_interfaces(coupling, x, if_buffer_send, if_buffer_recv);
}

double mpi_dotprod(index n, double* x, double* y) {
    double dotprod = 0;
    for (int i = 0; i < n; i++) {
        dotprod += x[i]*y[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotprod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return dotprod;
}