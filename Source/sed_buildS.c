#include "hpc.h"
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* compute local stiffness matrix
 * ordering w.r.t. [ p1, p2, p3] 
 * ax w.r.t. to i[3] = {0,0,1}, j[3] = {1,2,2}; */
void stima_laplace(double p1[2], double p2[2], double p3[2],  
                   index  typ, double dx[3], double ax[3])  
{
    int i, j;
    double d[3][2], fac;
    
    // Get vectors
    for (i = 0 ; i < 2 ; i++ ){
        d[0][i] = p3[i]-p2[i];
        d[1][i] = p1[i]-p3[i];
        d[2][i] = p2[i]-p1[i]; 
    }

    // Get factor for integral values (from integral transformation)
    fac = ( kappa(p1,typ) + kappa(p2,typ) + kappa(p3,typ) ) / 
          (6.0*(d[0][0]*d[1][1]-d[0][1]*d[1][0]));

    // Get super diagonal entries
    ax[0] = fac * (d[1][0]*d[0][0] + d[1][1]*d[0][1]); 
    ax[1] = fac * (d[2][0]*d[0][0] + d[2][1]*d[0][1]); 
    ax[2] = fac * (d[2][0]*d[1][0] + d[2][1]*d[1][1]); 

    // Get diagonal entries
    dx[0] = fac * (d[0][0]*d[0][0] + d[0][1]*d[0][1]); 
    dx[1] = fac * (d[1][0]*d[1][0] + d[1][1]*d[1][1]); 
    dx[2] = fac * (d[2][0]*d[2][0] + d[2][1]*d[2][1]); 
}


sed *sed_nz_pattern(mesh *M)                         
{
    index j, k, n, p, nC, nT,  nE, nz, *Elem, ind[6], *Si, *w, imin, imax;
    static int ai[3] = {0,0,1}, aj[3] = {1, 2, 2};

    // Declare matrix
    sed *S;
    
    // For convenience
    nT   = M->nelem; 
    nC   = M->ncoord; 
    Elem = M->elem; 

    // Get structure of sparse matrix
    n   = nC;
    nz  = n + 1 + 3 * nT;
    S   = sed_alloc (n, nz, 0) ;                         /* allocate result */ 
    w   = calloc (n, sizeof (index)) ;                   /* get workspace */
    if (!S || !w) return (sed_done (S, w, NULL, 0)) ;    /* out of memory */    
    Si = S->i ;  

    // Count entries per column
    for (k = 0 ; k < nT ; k++)                          
    {
        for (j = 0 ; j < 3 ; j++){
            ind[j] = Elem[7*k+j];                                   
        }
        for (j = 0 ; j < 3 ; j++)
        {
            w[ HPC_MIN( ind[ai[j]] , ind[aj[j]] ) ]++; /* off-diag. entries */
        }
    }

    // Column pointers
    hpc_cumsum(Si, w, n) ;

    // Offset
    for (k = 0 ; k < n ; k++){
        w[k] += n+1;                                 
    }

    // Offset column pointer
    for (k = 0 ; k < n+1 ; k++){
        Si[k] += n+1;                             
    }

    // Insert indeces
    for (k = 0 ; k < nT ; k++)                          
    {
        for (j = 0 ; j < 3 ; j++){
            ind[j] = Elem[7*k+j];                                   
        }
        for (j = 0 ; j < 3 ; j++){
            imin = HPC_MIN( ind[ai[j]] , ind[aj[j]] );
            imax = HPC_MAX( ind[ai[j]] , ind[aj[j]] );
            Si[w[imin]++] = imax ; 
        }
    }
    free(w);
    if (!sed_dupl(S)) return (sed_free(S)) ;          /* remove duplicates */
    return(S);
}


// Build the stiffness matrix in SED format
index sed_buildS(mesh *M, sed *T)                         
{
    index j, k, n, p, nC, nT, nz, *Elem, ind[3], *Ti, *w, imin, imax;
    
    // Indece aux arrays
    static int ai[3] = {0,0,1}, aj[3] = {1, 2, 2};

    // Element stiffness diagonal and super diagonal
    double dx[3], ax[3], *Coord, *Tx;
    
    nT      = M->nelem; 
    nC      = M->ncoord; 
    Coord   = M->coord; 
    Elem    = M->elem; 
    n       = T->n ; 
    Ti      = T->i ; 
    
    // Check if already allocated
    if (!(T->x)) Tx = T->x = calloc( Ti[n] , sizeof (double)) ;
    if (!Tx) return(0);
    
    // Get SED stiffness matrix entries
    for ( k = 0 ; k < nT; k++)
    {
        // Get node number per element
        for (j = 0 ; j < 3 ; j++){
            ind[j] = Elem[7*k+j];                                   
        }

        // Calculate local stiffness matrix for one element
        stima_laplace(Coord+2*ind[0],Coord+2*ind[1],Coord+2*ind[2],
                      Elem[7*k+7],dx,ax);

        // Add local stiffness diag to global
        for (j = 0 ; j < 3 ; j++){
            Tx[ind[j]] += dx[j];
        }

        // Fill-in column pointer and row indeces
        for (j = 0 ; j < 3 ; j++)
        {
            if (ind[ai[j]] < ind[aj[j]]){
                imin = ind[ai[j]]; 
                imax = ind[aj[j]]; 
            } else {
                imax = ind[ai[j]]; 
                imin = ind[aj[j]]; 
            }
    
            for (p = Ti[imin] ; p < Ti[imin+1] ; p++)
            {
                // Add local stiffness super diag to global
                if (Ti[p] == imax ) 
                {
                    Tx[p] += ax[j];
                    break;
                }
            }    
        }
    }
    return(1);
}
