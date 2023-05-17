#include "hpc.h"
#include <mpi.h>

coupling_data* mpi_split_interfaces(interface_data* interfaces, index* l2g, int n_nodes) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    coupling_data* coupling = malloc(sizeof(coupling_data));

    coupling->l2g = malloc(sizeof(index)*n_nodes);
    MPI_Datatype l2g_type;
    MPI_Type_contiguous(n_nodes, MPI_AINT, &l2g_type);
    MPI_Type_commit(&l2g_type);
    MPI_Scatter(l2g, 1, l2g_type, coupling->l2g, 1, l2g_type, 0, MPI_COMM_WORLD);

    // if (rank == 1) {
    //     printf("Hey I'm process %d and these are my l2g:\n", rank);
    //     for (int i = 0; i < n_nodes; i++) {
    //         printf("%zu\t", coupling->l2g[i]);
    //     }
    //     printf("\n");
    // }

    // Cross-Points
    coupling->nCrossPts = 4;
    index cross_pts[4] = {0, 1, 2, 3};
    coupling->crossPts = cross_pts;
    
    // Interfaces
    MPI_Scatter(rank == 0 ? interfaces->interface_counts : 0, 1, MPI_AINT, &coupling->ncoupl, 1, MPI_AINT, 0, MPI_COMM_WORLD);
    coupling->coupl = malloc(sizeof(index)*coupling->ncoupl*5);
    coupling->dcoupl = malloc(sizeof(index)*coupling->ncoupl);
    coupling->icoupl = malloc(sizeof(index*)*coupling->ncoupl);
    if (rank == 0) {
        index** interface_ids = interfaces->interface_ids;
        // copy information into local coupling data
        for (int j = 0; j < coupling->ncoupl; j++) {
            
            index interface_id = interfaces->interface_ids[0][j];
            coupling->dcoupl[j] = interfaces->dcoupl[interface_id];
            coupling->icoupl[j] = malloc(coupling->dcoupl[j]*sizeof(index));
            for (int k = 0; k < 5; k++) {
                coupling->coupl[5*j+k] = interfaces->coupl[5*interface_id+k];
            }
            for (int k = 0; k < coupling->dcoupl[j]; k++) {
                coupling->icoupl[j][k] = interfaces->icoupl[interface_id][k];
            }
        }
        // distribute information to other processes
        for (int i = 1; i < nof_p; i++) {
            for (int j = 0; j < interfaces->interface_counts[i]; j++) {
                index interface_id = interfaces->interface_ids[i][j];
                MPI_Send(&interfaces->coupl[5*interface_id], 5, MPI_AINT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&interfaces->dcoupl[interface_id], 1, MPI_AINT, i, 0, MPI_COMM_WORLD);
                MPI_Send(interfaces->icoupl[interface_id], interfaces->dcoupl[interface_id],
                    MPI_AINT, i, 0, MPI_COMM_WORLD);
            }
        }
    } else {
        for (int j = 0; j < coupling->ncoupl; j++) {
            MPI_Recv(&coupling->coupl[5*j], 5, MPI_AINT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&coupling->dcoupl[j], 1, MPI_AINT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            coupling->icoupl[j] = malloc(coupling->dcoupl[j]*sizeof(index));
            MPI_Recv(coupling->icoupl[j], coupling->dcoupl[j], MPI_AINT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    // if (rank == 0) {
    //     printf("Hey I'm process %d and this is my coupling data:\n", rank);
    //     printf("number interfaces: %zu\n", coupling->ncoupl);
    //     for (int i = 0; i < coupling->ncoupl; i++) {
    //         printf("Interface %d:\n", i);
    //         for (int j = 0; j < 5; j++) {
    //             printf("\t%zu", coupling->coupl[5*i+j]);
    //         }
    //         printf("\n");
    //         printf("\tnumber of coupling nodes: %zu\n", coupling->dcoupl[i]);
    //         for (int j = 0; j < coupling->dcoupl[i]; j++) {
    //             printf("\t%zu", coupling->icoupl[i][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    return coupling;
}