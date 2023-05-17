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

        // Save interfaces for debug
        char fname[200];
        sprintf(fname,"../Problem/rectangle_%dx%d",n_rows,n_cols);
        interface_data_write(interfaces, fname);
    }

    mesh* local_mesh = mesh_create_rect(1, 1, boundaries);
    local_mesh = mesh_multi_refine(local_mesh, refinements);
    int n_nodes = (int)local_mesh->ncoord;

    coupling_data* coupling = mpi_split_interfaces(interfaces, l2g_numbering, n_nodes);
    

    // TODO
    // speicher freigeben
    // colors
    // interfaces in lokalem numbering (und nicht globalem)
    
    MPI_Finalize();
}
