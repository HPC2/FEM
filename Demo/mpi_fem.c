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
  return ( 0.0 );
//  return ( x[0] * x[1] );
}

struct timeval tv[100];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))


char* jacobi = "jacobi";
char* gauss_seidel = "gs";
char* cg = "cg";
char* pcg = "pcg";

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    
    int n_rows;
    int n_cols;
    int refinements;
    index n_global_nodes;
    int n_global_crosspoints;
    char solver[10];
    char result_name[256];

    if (rank == 0) {
        if (argc != 6) {
          printf("Pass [rows] [cols] [refinements] [solver] [result_name] as arguments\n");
          return -1;
        }

        n_rows = atoi(argv[1]); //anzahl an Zeilen
        n_cols = atoi(argv[2]); //anzahl an Spalten
        refinements = atoi(argv[3]);
        strcpy(solver, argv[4]);
        strcpy(result_name, argv[5]);
        assert(nof_p == n_rows * n_cols);
    }

    TIME_SAVE(0);

    // broadcat number of refinements
    MPI_Bcast(&n_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&refinements, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(solver, 10, MPI_CHAR, 0, MPI_COMM_WORLD);

    index global_boundaries[4] = {1, 0, 0, 0};
    index* l2g_numbering;
    interface_data* interfaces;
    mesh* global_mesh;

    if (rank == 0) {
        // Global mesh
        global_mesh = mesh_create_rect(n_rows, n_cols, global_boundaries, 0.0, 0.0, 1.0, 1.0);
        
        if(!global_mesh) {
            printf("OOM\n");
            return 1;
        }

        // get number of global crosspoints (which is the number of nodes before refining)
        n_global_crosspoints = global_mesh->ncoord;

        // Save global mesh before refining for debug
        // char fname_glob_mesh_pre_ref[200];
        // sprintf(fname_glob_mesh_pre_ref,"../Problem/rectangle_%dx%d_global_0ref",n_rows,n_cols);
        // mesh_write(global_mesh, fname_glob_mesh_pre_ref);

        // Refine
        global_mesh = mesh_multi_refine(global_mesh, refinements);

        // get number of global nodes
        n_global_nodes = global_mesh->ncoord;

        // Save global mesh after refining for debug
        // char fname_glob_mesh_post_ref[200];
        // sprintf(fname_glob_mesh_post_ref,"../Problem/rectangle_%dx%d_global_%dref",n_rows,n_cols,refinements);
        // mesh_write(global_mesh, fname_glob_mesh_post_ref);

        // Numbering
        l2g_numbering = get_local_to_global_numbering(global_mesh, n_rows, n_cols, refinements);
        
        // write the local to global numbering
        //write_l2g(l2g_numbering, n_rows, n_cols, refinements);
        
        // Interfaces
        interfaces = rect_interface_data(global_mesh, n_rows, n_cols, refinements);

        // Save interfaces for debug
        // char fname_interfaces[200];
        // sprintf(fname_interfaces,"../Problem/rectangle_%dx%d",n_rows,n_cols);
        // interface_data_write(interfaces, fname_interfaces);
    }

    index* boundaries = mpi_boundaries(n_rows, n_cols, global_boundaries);

    double size_x = 1.0 / n_cols;
    double size_y = 1.0 / n_rows;
    double offset_x = (rank % n_cols) * size_x;
    double offset_y = (rank / n_cols) * size_y;

    mesh* local_mesh = mesh_create_rect(1, 1, boundaries, offset_x, offset_y, size_x, size_y);
    local_mesh = mesh_multi_refine(local_mesh, refinements);

    if (!strcmp(solver, gauss_seidel)) {
      mesh_flip_edge(local_mesh);
    }

    local_mesh->fixed = mesh_getFixed( local_mesh->ncoord, 
                                 local_mesh->bdry, 
                                 local_mesh->nbdry, 
                                &local_mesh->nfixed);

    index n    = local_mesh->ncoord;

    // broadcast number of global nodes
    MPI_Bcast(&n_global_nodes, 1, MPI_AINT, 0, MPI_COMM_WORLD);
    // create coupling data
    coupling_data* coupling = mpi_split_interfaces(interfaces, l2g_numbering, n, n_global_nodes);

    // coupling_data_print(coupling, rank);
    
    MPI_Bcast(&n_global_crosspoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
    coupling->n_global_cp = n_global_crosspoints;
    comm_buffers* buffers = alloc_comm_buffers(n_global_crosspoints, coupling);

    TIME_SAVE(1);

    // get pattern of matrix
    sed* A = sed_nz_pattern(local_mesh);
    if (!A) return(1);

    // Build stiffness matrix
    TIME_SAVE(10);
    if ( !sed_buildS(local_mesh, A) ) return(1); // assemble coefficient matrix
    MPI_Barrier(MPI_COMM_WORLD);
    TIME_SAVE(11);

    // Print the matrix
    // if (rank == 0) {
    //     sed_print(A, 1);
    // }

    // Get storage for rhs and solution
    // Initialize with zeros
    double* x    = calloc (n, sizeof(double));       // get workspace for sol
    double* b    = calloc (n, sizeof(double));       // get workspace for rhs

    // Build rhs (Volume and Neumann data)
    TIME_SAVE(20);
    mesh_buildRhs(local_mesh, b, F_vol, g_Neu); 
    MPI_Barrier(MPI_COMM_WORLD);
    TIME_SAVE(21);

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
    MPI_Barrier(MPI_COMM_WORLD);
    TIME_SAVE(22);

    index n_iter;

    TIME_SAVE(30);
    if (!strcmp(solver, jacobi)) {
      n_iter = mpi_jacobi(A, coupling, buffers, local_mesh, x, b);
    } else if (!strcmp(solver, gauss_seidel)) {
      mpi_gs(A, coupling, buffers, local_mesh, x, b);
    } else if (!strcmp(solver, cg)) {
      n_iter = mpi_cg(A, coupling, buffers, local_mesh, x, b);
    } else if (!strcmp(solver, pcg)) {
      n_iter = mpi_pcg(A, coupling, buffers, local_mesh, x, b);
    } else {
      printf("Solver unknown\n");
      MPI_Finalize();
      return -1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    TIME_SAVE(31);

    // Write result
    index n_iter_max=0;
    index n_iter_min=0;
    index n_iter_sum=0;
    MPI_Reduce(&n_iter, &n_iter_max,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&n_iter, &n_iter_min,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(&n_iter, &n_iter_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    if (rank == 0) {
      result_write(
        result_name,
        n_rows*n_cols,
        n_rows,
        n_cols,
        refinements,
        solver,
        n_iter_max,
        n_iter_min,
        n_iter_sum,
        n_global_nodes,
        global_boundaries,
        (int)TIME_ELAPSED(0, 1),
        (int)TIME_ELAPSED(10,11),
        (int)TIME_ELAPSED(20,22),
        (int)TIME_ELAPSED(30,31)
      );
    }


    // save A, x and b 
    //double* global_x = mpi_assemble_t1_vec(coupling, x, n);
    // double* global_rhs = mpi_assemble_t2_vec(coupling, b, n);
    // double* global_A   = mpi_assemble_A(A, coupling);
    if (rank == 0) {
        //char buf[200];
        //sprintf(buf, "../Problem/x_mpi_%s", solver);
        //print_dmatrix(global_x, n_global_nodes, 1, false, buf, "dat");
        // print_dmatrix(global_rhs, n_global_nodes, 1, false, "../Problem/b-test", "dat");
        // print_dmatrix(global_A, n_global_nodes, n_global_nodes, true, "../Problem/A-test", "dat");
    }
    
    // free some stuff
    mesh_free(local_mesh);
    sed_free(A);
    if (rank == 0) mesh_free(global_mesh);
    free(x);
    free(b); 
    free_comm_buffers(buffers);
    
    // stop MPI and return
    MPI_Finalize();
    return 0;
}
