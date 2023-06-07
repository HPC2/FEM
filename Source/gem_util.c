#include "hpc.h"

// Allocate a gem matrix
gem *gem_alloc(index m, index n, index incRow, index incCol)
{
    gem *A = (gem*) malloc(sizeof (gem)) ;    /* allocate the gem struct */
    if (!A) {
        printf("OOM!\n");
        return (NULL) ;                   /* out of memory */
    } 

    // Assign the input values
    A->n      = n ;                         /* define dimensions */
    A->m      = m ;                         /* define dimensions */
    A->incRow = incRow ;                         /* row increment */
    A->incCol = incCol ;                         /* col increment*/

    // Allocate storage for entries
    A->x      = (double*) calloc(n*m, sizeof(double)) ;
    if (!A->x) {
        printf("OOM!\n");
        return (NULL);
    }

    return A;
}

// free a dense matrix 
gem *gem_free(gem *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->x) ;
    free (A);
    return (NULL) ; 
}
