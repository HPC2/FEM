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
//  return ( 0.0 );
  return ( x[0] * x[1] );
}


int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Pass [rows] [cols] [refinements] as arguments\n");
        return -1;
    }

    int n_rows = atoi(argv[1]); //anzahl an Zeilen
    int n_cols = atoi(argv[2]); //anzahl an Spalten
    int refinements = atoi(argv[3]);

    char pdir[] = "../Problem/";

    index boundaries[4] = {0, 0, 0, 0};
    mesh* Mesh = mesh_create_rect(n_rows, n_cols, boundaries, 0.0, 0.0);
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

    // mesh_print(Mesh, 0);
    // char fname[64]; sprintf(fname, "%srectangle_%dx%d_refined%dtimes", pdir, n_rows, n_cols, refinements);
    // mesh_write(Mesh, fname);

    // get pattern of matrix
    sed* A = sed_nz_pattern(Mesh);
    if (!A) return(1);

    // Build stiffness matrix
    if ( !sed_buildS(Mesh, A) ) return(1); // assemble coefficient matrix

    // Get storage for rhs and solution
    index n    = A->n;
    double* b    = calloc (n, sizeof(double));       // get workspace for rhs

    // Build rhs (Volume and Neumann data)
    mesh_buildRhs(Mesh, b, F_vol, g_Neu); 

    gem* A_full = sed_to_dense(A, true);
    print_dmatrix(b, n, 1, false, "../Problem/rhs", "dat");
    // print_dmatrix(A_full->x, n, n, true, "../Problem/A", "dat");
    // printf("global rhs:\n");
    // for (index i = 0; i < Mesh->ncoord; i++) {
    //   printf("%4.4lf\n", b[i]);
    // }

    mesh_free(Mesh);
}