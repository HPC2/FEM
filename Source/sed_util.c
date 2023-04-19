#include "hpc.h"
// allocate a sparse square matrix with extracted diagonal in 
// compressed col. form
sed *sed_alloc(index n, index nzmax, index values)
{
    sed *A = calloc (1, sizeof (sed)) ;    /* allocate the sed struct */
    if (!A) return (NULL) ;                /* out of memory */

    A->n     = n ;                         /* define dimensions and nzmax */
    A->nzmax = nzmax = HPC_MAX (nzmax, n+1) ;

    A->i     = malloc (nzmax * sizeof (index)) ;
    if (!A->i) return (NULL);

    A->i[A->n] = nzmax;
    A->x = values ? malloc(nzmax * sizeof (double)) : NULL ;
    return ((!A->i || (values && !A->x)) ? sed_free (A) : A) ;
}

/* change the max # of entries sparse matrix */
index sed_realloc(sed *A, index nzmax)
{
    index ok, oki, okx = 1 ;
    if (!A) return 0 ;
    if (nzmax <= 0){
        nzmax = A->i[A->n];
    }
    nzmax = HPC_MAX(nzmax, 1+A->n) ;
    A->i  = hpc_realloc(A->i, nzmax, sizeof (index), &oki) ;
    if (A->x){
        A->x = hpc_realloc(A->x, nzmax, sizeof (double), &okx) ;
    }

    ok = (oki && okx) ;
    if (ok){
        A->nzmax = nzmax ;
    }
    return ok ;
}

/* free a sparse matrix */
sed *sed_free(sed *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->i) ;      /* free the sed struct and return NULL */
    free (A->x) ;
    free (A);
    return (NULL) ; 
}

/* free workspace and return a sparse matrix result */
sed *sed_done(sed *C, void *w, void *x, index ok)
{
    free (w) ;                       /* free workspace */
    free (x) ;
    return (ok ? C : sed_free (C)) ;   /* return result if OK, else free it */
}


// sed matrix to dense
gem *sed_to_dense(const sed *A, bool sym)
{
    index m, n, incRowB, incColB;

    // Assign dimension
    m = n = A->n;

    // Col major
    incRowB = 1;
    incColB = n;

    // Get storage for dense B
    gem* B = gem_alloc(m, n, incRowB, incColB);

    // Diagonal
    for (size_t i = 0; i < n; ++i){
        B->x[i*incRowB + i*incColB] = A->x[i];
    }

    // Super-/Sub-Diagonal
    if (sym){
        for (size_t j=0; j < n; ++j){
            for (size_t p=A->i[j]; p < A->i[j+1]; ++p){
                B->x[A->i[p]*incRowB + j*incColB] = A->x[p];
                B->x[j*incRowB + A->i[p]*incColB] = A->x[p];
            }
        }
    } else {
        for (size_t j=0; j < n; ++j){
            for (size_t p=A->i[j]; p < A->i[j+1]; ++p){
                B->x[A->i[p]*incRowB + j*incColB] = A->x[p];
            }
        }
    }

    return B;
}





