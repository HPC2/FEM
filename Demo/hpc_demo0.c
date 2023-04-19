#include "hpc.h"
#include <inttypes.h>
#include <time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

double kappa( double x[2], index typ )
{
  return ( 1.0 );
}

double F_vol( double x[2], index typ )
{
  return ( 0.0 );
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


int main (int argc, char **argv)
{
    // Number of refinements specified by user input
    index N = 0;
    // Relative path of initial mesh
    char *Pdir = "../Problem/";
    // Problem name specified by user, e.g. problem1 
    char fname[64];

    // Declare some aux variables
    index n, k, ncoord, nelem, nbdry, nfixed, nedges, *fixed, total, *bdry,
          cnt = 0, max_iter = 5000,
          brief_print = 1;

    double *b, *x, *w, *resi, *Coord, x1[2], resi_norm = 0, max_tol=1e-8;
 
    // Exception handling
    printf("\n========================================\n");
    if (argc < 2 )
    { 
        printf("Problem not specified\n"); 
        return(1); 
    } 

    // Get problem as parameter from user
    sprintf(fname,"%s%s",Pdir,argv[1]); 

    // Get no. of refinements from user
    if (argc>2)
    {                        
        // Min and max number of refinements
        if ( (atoi(argv[2]) > 0) & (atoi(argv[2]) < 13))
        {
            N = atoi(argv[2]);
        }
    } 

    // Declare mesh and sed-matrix pointer
    mesh *H[2];
    sed  *A   ;

    // Print basic information
    printf("Load data form %s, no. refinements = %g\n", fname, (double) N);   

   
    TIME_SAVE(0);
    // Load specific problem
    // Load geometry
    H[0] = mesh_load(fname);

    // Get the edge to node information (necessary for refine)
    mesh_getEdge2no( H[0]->nelem, 
                     H[0]->elem, 
                    &H[0]->nedges,
                    &H[0]->edge2no);

    // Print mesh information
    mesh_print(H[0], brief_print); 

    // Refine mesh N times
    //printf("======== Start %3g refinements ========\n",(double)  N);
    for (size_t i=0; i<N; ++i)
    {
        // Refine
        H[1] = mesh_refine(H[0]);

        // Get the edge to node information (necessary for refine)
        mesh_getEdge2no( H[1]->nelem, 
                         H[1]->elem, 
                        &H[1]->nedges,
                        &H[1]->edge2no);

        if (mesh_free(H[0]) != NULL) {return(0);}
        H[0] = H[1];
    }

    // Get the fixed node information
    H[0]->fixed = mesh_getFixed( H[0]->ncoord, 
                                 H[0]->bdry, 
                                 H[0]->nbdry, 
                                &H[0]->nfixed);


    // Temp mesh pointer is not needed anymore
    H[1] = NULL;

    // Print refined mesh information
    mesh_print(H[0], brief_print); 

    // get pattern of matrix
    A = sed_nz_pattern(H[0]);
    if (!A) return(1);

    // Build stiffness matrix
    if ( !sed_buildS(H[0],A) ) return(1); // assemble coefficient matrix
    TIME_SAVE(1);

    // Print the matrix
    sed_print(A, brief_print);

    // Get storage for rhs and solution
    n    = A->n;
    // Initialize with zeros
    x    = calloc (n, sizeof(double));       // get workspace for sol
    w    = calloc (n, sizeof(double));       // get temporary workspace
    b    = calloc (n, sizeof(double));       // get workspace for rhs
    resi = calloc (n, sizeof(double));       // get workspace for residual

    // Build rhs (Volume and Neumann data)
    mesh_buildRhs(H[0], b, F_vol, g_Neu); 
    TIME_SAVE(2);

    // For convenience
    ncoord      = H[0]->ncoord ; 
    nelem       = H[0]->nelem ; 
    nbdry       = H[0]->nbdry ; 
    nfixed      = H[0]->nfixed ; 
    nedges      = H[0]->nedges; 

    fixed       = H[0]->fixed ; 
    bdry        = H[0]->bdry; 
    Coord       = H[0]->coord;

    // Adjust rhs to incorporate Dirichlet data
    // x <-- u_D (at dirichlet nodes)
    for ( k = 0; k < nfixed; ++k)
    {
        x1[0] = Coord[2 * fixed[k]]; 
        x1[1] = Coord[2 * fixed[k]+1];

        x[fixed[k]] = u_D(x1);
    }

    // Solve with sym. Gauss-Seidel iterations (and don't touch dirichlet nodes)
    for (k = 0; k<max_iter; ++k){

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
        resi_norm = ddot(A->n, resi, 1, resi, 1);
        if (resi_norm < max_tol){
            break;
        }

        // Sym Gauss-Seidel iterations
        sed_gs_constr(A, b, x, w, H[0]->fixed, H[0]->nfixed, 1); 
        sed_gs_constr(A, b, x, w, H[0]->fixed, H[0]->nfixed, 0); 
    }

    TIME_SAVE(3);

    // print times
    printf("\n");
    printf("Time load & refine mesh = %9i ns\n", (int) TIME_ELAPSED(0,1));
    printf("Time building rhs            = %9i ns\n", (int) TIME_ELAPSED(1,2));
    printf("Time Solve                   = %9i ns\n", (int) TIME_ELAPSED(2,3));

    // print residual
    dcopy(A->n, b, 1, resi, 1);
    sysed_spmv(-1, A, x, 1, 1, resi, 1);
    for ( size_t i = 0; i < nfixed; ++i){
        resi[fixed[i]] = 0;
    }
    resi_norm = ddot(A->n, resi, 1, resi, 1);
    printf("\n|| A*x_k - b ||_2^2 = %10g\n", resi_norm);

    printf("\n========================================\n");

    // Free everything
    mesh_free(H[1]); 
    sed_free(A);
    free(x);
    free(b);
    free(w);
    free(resi);

}


