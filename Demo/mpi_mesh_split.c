#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "hpc.h"

double kappa( double x[2], index typ )
{
  return ( 1.0 );
}

double F_vol( double x[2], index typ )
{
  return ( 1.0 );
}

double g_Neu( double x[2], index typ )
{
  return ( x[0] * x[1] );
}

double u_D( double x[2])
{
//  return ( 0.0 );
  return ( x[0] * x[1] );
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    
    int n_rows;
    int n_cols;
    int refinements;
    index n_global_nodes;
    int n_global_crosspoints;

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

    // broadcat number of refinements
    MPI_Bcast(&n_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&refinements, 1, MPI_INT, 0, MPI_COMM_WORLD);

    index global_boundaries[4] = {0, 0, 0, 0};
    index* l2g_numbering;
    interface_data* interfaces;
    mesh* global_mesh;

    if (rank == 0) {
        // Global mesh
        global_mesh = mesh_create_rect(n_rows, n_cols, global_boundaries, 0.0, 0.0);
        
        if(!global_mesh) {
            printf("OOM\n");
            return 1;
        }

        // get number of global crosspoints (which is the number of nodes before refining)
        n_global_crosspoints = global_mesh->ncoord;

        // Save global mesh before refining for debug
        char fname_glob_mesh_pre_ref[200];
        sprintf(fname_glob_mesh_pre_ref,"../Problem/rectangle_%dx%d_global_0ref",n_rows,n_cols);
        mesh_write(global_mesh, fname_glob_mesh_pre_ref);

        // Refine
        global_mesh = mesh_multi_refine(global_mesh, refinements);

        // get number of global nodes
        n_global_nodes = global_mesh->ncoord;

        // Save global mesh after refining for debug
        char fname_glob_mesh_post_ref[200];
        sprintf(fname_glob_mesh_post_ref,"../Problem/rectangle_%dx%d_global_%dref",n_rows,n_cols,refinements);
        mesh_write(global_mesh, fname_glob_mesh_post_ref);

        // Numbering
        l2g_numbering = get_local_to_global_numbering(global_mesh, n_rows, n_cols, refinements);
        
        // write the local to global numbering
        write_l2g(l2g_numbering, n_rows, n_cols, refinements);
        
        // Interfaces
        interfaces = rect_interface_data(global_mesh, n_rows, n_cols, refinements);

        // Save interfaces for debug
        char fname_interfaces[200];
        sprintf(fname_interfaces,"../Problem/rectangle_%dx%d",n_rows,n_cols);
        interface_data_write(interfaces, fname_interfaces);
    }


    index* boundaries = mpi_boundaries(n_rows, n_cols, global_boundaries);

    double offset_x = rank % n_cols;
    double offset_y = rank / n_cols;

    mesh* local_mesh = mesh_create_rect(1, 1, boundaries, offset_x, offset_y);
    local_mesh = mesh_multi_refine(local_mesh, refinements);

    local_mesh->fixed = mesh_getFixed( local_mesh->ncoord, 
                                 local_mesh->bdry, 
                                 local_mesh->nbdry, 
                                &local_mesh->nfixed);

    // get pattern of matrix
    sed* A = sed_nz_pattern(local_mesh);
    if (!A) return(1);

    // Build stiffness matrix
    if ( !sed_buildS(local_mesh, A) ) return(1); // assemble coefficient matrix

    // Print the matrix
    if (rank == 0) {
        // sed_print(A, 1);
    }

    // Get storage for rhs and solution
    index n    = A->n;
    // Initialize with zeros
    double* x    = calloc (n, sizeof(double));       // get workspace for sol
    double* b    = calloc (n, sizeof(double));       // get workspace for rhs

    // Build rhs (Volume and Neumann data)
    mesh_buildRhs(local_mesh, b, F_vol, g_Neu); 

    // broadcast number of global nodes
    MPI_Bcast(&n_global_nodes, 1, MPI_AINT, 0, MPI_COMM_WORLD);
    // create coupling data
    coupling_data* coupling = mpi_split_interfaces(interfaces, l2g_numbering, n, n_global_nodes);

    // coupling_data_print(coupling, rank);
    
    MPI_Bcast(&n_global_crosspoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
    coupling->n_global_cp = n_global_crosspoints;
    comm_buffers* buffers = alloc_comm_buffers(n_global_crosspoints, coupling);

    index* fixed = local_mesh->fixed;
    index nfixed = local_mesh->nfixed;
    double* Coord = local_mesh->coord;
    double x1[2];
    for ( int k = 0; k < nfixed; ++k)
    {
        x1[0] = Coord[2 * fixed[k]]; 
        x1[1] = Coord[2 * fixed[k]+1];

        x[fixed[k]] = u_D(x1);
    }

    mpi_jacobi(A, coupling, buffers, local_mesh, x, b);

    // save x
    double* global_x = mpi_assemble_t1_vec(coupling, x, n);
    if (rank = 0) {
        print_dmatrix(global_x, n_global_nodes, 1, false, "../Problem/x-test1", "dat");
    }

    if (rank == 0) {
        mesh_free(global_mesh);
    }
    MPI_Finalize();
    return 0;
}
