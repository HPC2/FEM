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

    Mesh->fixed = mesh_getFixed( Mesh->ncoord, 
                                 Mesh->bdry, 
                                 Mesh->nbdry, 
                                &Mesh->nfixed);

    // get pattern of matrix
    sed* A = sed_nz_pattern(Mesh);
    if (!A) return(1);

    // Build stiffness matrix
    if ( !sed_buildS(Mesh, A) ) return(1); // assemble coefficient matrix

    // Get storage for rhs and solution
    index n    = A->n;
    // Initialize with zeros
    double* x    = calloc (n, sizeof(double));       // get workspace for sol
    double* w    = calloc (n, sizeof(double));       // get temporary workspace
    double* b    = calloc (n, sizeof(double));       // get workspace for rhs
    double* resi = calloc (n, sizeof(double));       // get workspace for residual

    // Build rhs (Volume and Neumann data)
    mesh_buildRhs(Mesh, b, F_vol, g_Neu);

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

    // Solve with sym. Gauss-Seidel iterations (and don't touch dirichlet nodes)
    for (index k = 0; k< 1000; ++k){

        // Calculate square norm of residual
        // resi <-- b
        dcopy(A->n, b, 1, resi, 1);
        // resi <-- b - A*x
        sysed_spmv(-1, A, x, 1, 1, resi, 1);
        // set fixed nodes to zero
        for ( size_t i = 0; i < nfixed; ++i){
            resi[fixed[i]] = 0;
        }
        // resi_norm <-- || resi ||_2^2
        double resi_norm = ddot(A->n, resi, 1, resi, 1);
        if (resi_norm < 1e-5){
            break;
        }

        // Sym Gauss-Seidel iterations
        sed_gs_constr(A, b, x, w, fixed, nfixed, 1); 
        sed_gs_constr(A, b, x, w, fixed, nfixed, 0); 
    }

    print_dmatrix(x, n, 1, false, "../Problem/x", "dat");
    print_dmatrix(b, n, 1, false, "../Problem/b", "dat"); 
    gem* A_full = sed_to_dense(A, true);
    double* A_data = A_full->x;
    print_dmatrix(A_data, n,n, true, "../Problem/A", "dat");

    mesh_free(Mesh);
}