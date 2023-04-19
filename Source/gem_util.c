#include "hpc.h"

// Allocate a gem matrix
gem *gem_alloc(index m, index n, index incRow, index incCol)
{
    gem *A = (gem*) malloc(sizeof (gem)) ;    /* allocate the gem struct */
    if (!A) return (NULL) ;                   /* out of memory */

    // Assign the input values
    A->n      = n ;                         /* define dimensions */
    A->m      = m ;                         /* define dimensions */
    A->incRow = incRow ;                         /* row increment */
    A->incCol = incCol ;                         /* col increment*/

    // Allocate storage for entries
    A->x      = (double*) malloc(n*m*sizeof(double)) ;
    if (!A->x) return (NULL);

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


