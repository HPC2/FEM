#include "hpc.h"
/* print a sparse matrix; use %g for integers to avoid differences with index */
index sed_print (const sed *A, index brief)
{
    index p, j, n, nzmax, *Ai ;
    double *Ax ;

    // Exception handling
    if (!A) { printf ("(null)\n") ; return (0) ; }

    // For convenience
    n       = A->n ; 
    Ai      = A->i ; 
    Ax      = A->x ; 
    nzmax   = A->nzmax ;
    
    // Print basic information
    printf ("%g-by-%g, nzmax: %g nnz: %g\n", (double) n,
            (double) n, (double) nzmax, (double) Ai[n]) ;

    // Print diagonal
    printf ("Diagonal entries:\n"); 
    for (j = 0 ; j < n ; j++)
    {
        printf ("      %g : %g\n", (double) j, Ax ? Ax[j] : 1) ;
        if (brief && j > 10) { 
            printf ("  ...\n") ; 
            break ; 
        }
    }

    // Print off-diagonal
    printf ("Off-diagonal entries:\n"); 
    for (j = 0 ; j < n ; j++)
    {
        printf("    col %g : locations %g to %g\n", 
               (double) j, 
               (double) (Ai[j]), 
               (double) (Ai[j+1]-1)) ;

        for (p = Ai[j] ; p < Ai[j+1] ; p++)
        {
            printf("(%3g,%3g) : %g\n", 
                   (double) (Ai[p]), 
                   (double) j, 
                   Ax[p]) ;

            if (brief && j > 10) { 
                printf("  ...\n") ; 
                return (1) ; 
            }
        }
    }

    return (1) ;
}

