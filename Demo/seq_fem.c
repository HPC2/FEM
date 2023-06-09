#include <stdio.h>
#include <stdlib.h>
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
//   return ( x[0] * x[1] );
}

struct timeval tv[100];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

char* jacobi = "jacobi";
char* gauss_seidel = "gs";
char* cg = "cg";
char* pcg = "pcg";


int main(int argc, char **argv) {
    if (argc != 6) {
        printf("Pass [rows] [cols] [refinements] [solver] [result_name] as arguments\n");
        return -1;
    }

    int n_rows = atoi(argv[1]); //anzahl an Zeilen
    int n_cols = atoi(argv[2]); //anzahl an Spalten
    int refinements = atoi(argv[3]);
    char* solver = argv[4];
    char* result_name = argv[5];

    TIME_SAVE(0);
    index boundaries[4] = {1, 0, 0, 0};
    mesh* Mesh = mesh_create_rect(n_rows, n_cols, boundaries, 0.0, 0.0, 1.0, 1.0);
    if(!Mesh) {
        printf("OOM\n");
        return 1;
    }

    mesh_getEdge2no( Mesh->nelem, 
                     Mesh->elem, 
                    &Mesh->nedges,
                    &Mesh->edge2no);
    

    // refine the mesh
    for(int i=0; i<refinements; i++) {
        Mesh = mesh_refine(Mesh);
        // update edge/node info??
        mesh_getEdge2no( Mesh->nelem, 
                     Mesh->elem, 
                    &Mesh->nedges,
                    &Mesh->edge2no);
    }

    // mesh_flip_edge(Mesh);

    Mesh->fixed = mesh_getFixed( Mesh->ncoord, 
                                 Mesh->bdry, 
                                 Mesh->nbdry, 
                                &Mesh->nfixed);
    TIME_SAVE(1);
    
    // get pattern of matrix
    sed* A = sed_nz_pattern(Mesh);
    if (!A) return(1);

    // Build stiffness matrix
    TIME_SAVE(10);
    if ( !sed_buildS(Mesh, A) ) return(1); // assemble coefficient matrix
    TIME_SAVE(11);

    // Get storage for rhs and solution
    index n    = A->n;
    // Initialize with zeros
    double* x    = calloc (n, sizeof(double));       // get workspace for sol
    double* b    = calloc (n, sizeof(double));       // get workspace for rhs

    // Build rhs (Volume and Neumann data)
    TIME_SAVE(20);
    mesh_buildRhs(Mesh, b, F_vol, g_Neu);
    TIME_SAVE(21);

    // For convenience
    index nfixed = Mesh->nfixed ; 

    index* fixed       = Mesh->fixed ; 
    double* Coord       = Mesh->coord;
    double x1[2];
    // Adjust rhs to incorporate Dirichlet data
    // x <-- u_D (at dirichlet nodes)
    for ( index k = 0; k < nfixed; ++k)
    {
        x1[0] = Coord[2 * fixed[k]]; 
        x1[1] = Coord[2 * fixed[k]+1];

        x[fixed[k]] = u_D(x1);
    }
    TIME_SAVE(22);

    index n_iter;

    TIME_SAVE(30);
    if (!strcmp(solver, jacobi)) {
      n_iter = seq_jacobi(A, Mesh, x, b);
    } else if (!strcmp(solver, gauss_seidel)) {
      n_iter = seq_gs(A, Mesh, x, b);
    } else if (!strcmp(solver, cg)) {
      n_iter = seq_cg(A, Mesh, x, b);
    } else if (!strcmp(solver, pcg)) {
      n_iter = seq_pcg(A, Mesh, x, b);
    } else {
      printf("Solver unknown\n");
      return -1;
    }
    TIME_SAVE(31);

    result_write(
      result_name,
      1,
      n_rows,
      n_cols,
      refinements,
      solver,
      n_iter,
      n_iter,
      n_iter,
      Mesh->ncoord,
      boundaries,
      (int)TIME_ELAPSED(0, 1),
      (int)TIME_ELAPSED(10,11),
      (int)TIME_ELAPSED(20,22),
      (int)TIME_ELAPSED(30,31),
      0,
      0
    );

    //char buf[200];
    //sprintf(buf, "../Problem/x_seq_%s", solver);

    //print_dmatrix(x, n, 1, false, buf, "dat");
    
    // print_dmatrix(b, n, 1, false, "../Problem/b", "dat"); 
    // gem* A_full = sed_to_dense(A, true);
    // double* A_data = A_full->x;
    // print_dmatrix(A_data, n,n, true, "../Problem/A", "dat");

    mesh_free(Mesh);
    free(x);
    free(b);
    sed_free(A);
}