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

    // if (rank == 0) {
    //     printf("Rank %d boundaries: %td, %td, %td, %td\n", rank, boundaries[0], boundaries[1], boundaries[2], boundaries[3]);
    // }
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
    double* w    = calloc (n, sizeof(double));       // get temporary workspace
    double* b    = calloc (n, sizeof(double));       // get workspace for rhs
    double* resi = calloc (n, sizeof(double));       // get workspace for residual

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

    // jacobi solver
    double* inv_diag = accumulate_inv_diag(coupling, A, buffers);
    dcopy(n, b, 1, resi, 1); // resi <- b
    sysed_spmv(-1, A, x, 1, 1, resi, 1);  // resi <- resi - A*x
    dcopy(n, resi, 1, w, 1); // w <- resi
    mpi_convert_type2_to_type1(coupling, w, buffers);  // convert w
    dmult(n, inv_diag, 1, w, 1); // w <- D^-1 * w
    double sigma = mpi_dotprod(n, w, resi);
    printf("rank %d: sigma = %lf\n", rank, sigma);
    double sigma0 = sigma;
    double tol = 1e-6;
    double omega = 1; // full jacobi steps for now
    while (sigma > tol*sigma0) {
      // update x
      daxpy(n, omega, w, 1, x, 1);
      // recalc residuum
      dcopy(n, b, 1, resi, 1); // resi <- b
      sysed_spmv(-1, A, x, 1, 1, resi, 1);  // resi <- resi - A*x
      for ( size_t i = 0; i < nfixed; ++i){
            resi[fixed[i]] = 0;
      }
      // recalc w
      dcopy(n, resi, 1, w, 1); // w <- resi
      mpi_convert_type2_to_type1(coupling, w, buffers);  // convert w
      dmult(n, inv_diag, 1, w, 1); // w <- D^-1 * w
      // recalc error 
      sigma = mpi_dotprod(n, w, resi);
      if (rank == 0) {
        printf("sigma = %lf\n", sigma);
        for (int i = 0; i < n; i++) {
          printf("%lf\t", x[i]);
        }
        printf("\n");
      }
    }

    // save x
    double* global_x = mpi_assemble_t1_vec(coupling, x, n);
    if (rank = 0) {
        print_dmatrix(global_x, n_global_nodes, 1, false, "../Problem/x-test1", "dat");
    }

    // double* global_b2 = mpi_assemble_t2_vec(coupling, b, n_nodes); // warning: global_b only contains smth sensible for rank 0
    // mpi_convert_type2_to_type1(coupling, b, cp_buffer, if_buffer1, if_buffer2);
    // double* global_b1 = mpi_assemble_t1_vec(coupling, b, n_nodes);
    // if (rank == 0) {
    //   print_dmatrix(global_b2, n_global_nodes, 1, false, "../Problem/rhs-test2", "dat");
    //   print_dmatrix(global_b1, n_global_nodes, 1, false, "../Problem/rhs-test1", "dat");
    //   for (int i = 0; i < n_global_nodes; i++) {
    //     assert(fabs(global_b1[i] - global_b2[i]) < 1e-13);
    //   }
    // }

    if (rank == 0) {
        // mesh_free(global_mesh);
    }
    MPI_Finalize();
    return 0;
}
