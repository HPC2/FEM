#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "hpc.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    
    int n_rows;
    int n_cols;
    int refinements;

    if (rank == 0) {
        if (argc != 4) {
            printf("Pass [rows] [cols] [refinments] as arguments\n");
            return -1;
        }

        n_rows = atoi(argv[1]); //anzahl an Zeilen
        n_cols = atoi(argv[2]); //anzahl an Spalten
        refinements = atoi(argv[3]);
        assert(nof_p == n_rows * n_cols);
    }

    MPI_Bcast(&refinements, 1, MPI_INT, 0, MPI_COMM_WORLD);
    bool boundaries[4] = {0, 0, 0, 0};
    index* l2g_numbering;
    interface_data* interfaces;

    if (rank == 0) {
        
        mesh* global_mesh = mesh_create_rect(n_rows, n_cols, boundaries);
        if(!global_mesh) {
            printf("OOM\n");
            return 1;
        }

        global_mesh = mesh_multi_refine(global_mesh, refinements);

        interfaces = rect_interface_data(global_mesh, n_rows, n_cols, refinements);
        l2g_numbering = get_local_to_global_numbering(global_mesh, n_rows, n_cols, refinements);
        mesh_free(global_mesh);
    }

    mesh* local_mesh = mesh_create_rect(1, 1, boundaries);
    local_mesh = mesh_multi_refine(local_mesh, refinements);
    int n_nodes = (int)local_mesh->ncoord;

    coupling_data coupling;
    // L2G
    coupling.l2g = malloc(sizeof(index)*n_nodes);
    MPI_Datatype l2g_type;
    MPI_Type_contiguous(n_nodes, MPI_AINT, &l2g_type);
    MPI_Type_commit(&l2g_type);
    MPI_Scatter(l2g_numbering, 1, l2g_type, coupling.l2g, 1, l2g_type, 0, MPI_COMM_WORLD);

    if (rank == 1) {
        printf("Hey I'm process %d and these are my l2g:\n", rank);
        for (int i = 0; i < n_nodes; i++) {
            printf("%zu\t", coupling.l2g[i]);
        }
        printf("\n");
    }

    // Cross-Points
    coupling.nCrossPts = 4;
    index cross_pts[4] = {0, 1, 2, 3};
    coupling.crossPts = cross_pts;
    
    // Interfaces
    MPI_Scatter(rank == 0 ? interfaces->interface_counts : 0, 1, MPI_AINT, &coupling.ncoupl, 1, MPI_AINT, 0, MPI_COMM_WORLD);
    coupling.coupl = malloc(sizeof(index)*coupling.ncoupl*5);
    coupling.dcoupl = malloc(sizeof(index)*coupling.ncoupl);
    coupling.icoupl = malloc(sizeof(index*)*coupling.ncoupl);
    if (rank == 0) {
        index** interface_ids = interfaces->interface_ids;
        // copy information into local coupling data
        for (int j = 0; j < coupling.ncoupl; j++) {
            
            index interface_id = interfaces->interface_ids[0][j];
            coupling.dcoupl[j] = interfaces->dcoupl[interface_id];
            coupling.icoupl[j] = malloc(coupling.dcoupl[j]*sizeof(index));
            for (int k = 0; k < 5; k++) {
                coupling.coupl[5*j+k] = interfaces->coupl[5*interface_id+k];
            }
            for (int k = 0; k < coupling.dcoupl[j]; k++) {
                coupling.icoupl[j][k] = interfaces->icoupl[interface_id][k];
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
        for (int j = 0; j < coupling.ncoupl; j++) {
            MPI_Recv(&coupling.coupl[5*j], 5, MPI_AINT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&coupling.dcoupl[j], 1, MPI_AINT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            coupling.icoupl[j] = malloc(coupling.dcoupl[j]*sizeof(index));
            MPI_Recv(coupling.icoupl[j], coupling.dcoupl[j], MPI_AINT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    if (rank == 0) {
        printf("Hey I'm process %d and this is my coupling data:\n", rank);
        printf("number interfaces: %zu\n", coupling.ncoupl);
        for (int i = 0; i < coupling.ncoupl; i++) {
            printf("Interface %d:\n", i);
            for (int j = 0; j < 5; j++) {
                printf("\t%zu", coupling.coupl[5*i+j]);
            }
            printf("\n");
            printf("\tnumber of coupling nodes: %zu\n", coupling.dcoupl[i]);
            for (int j = 0; j < coupling.dcoupl[i]; j++) {
                printf("\t%zu", coupling.icoupl[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    // TODO
    // speicher freigeben
    // colors
    // interfaces in lokalem numbering (und nicht globalem)
    
    MPI_Finalize();
}
