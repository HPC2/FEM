#include "hpc.h"

void set_dir(const mesh *M, double *u ){
    
    index j, k, *Bdry, nB, ind[2];
    
    nB      = M->nbdry; 
    Bdry    = M->bdry;

    for ( k = 0 ; k < nB; k++)
    {
        if (Bdry[4*k+3] == 0) 
        {
            ind[0] = Bdry[4*k+0];
            ind[1] = Bdry[4*k+1];

            for (j = 0 ; j < 2 ; j++){
                u[ind[j]] = 0;  
            }
        }
    } 


}

void set_dir_u_D(const mesh *M, double *u , const double* u_D ){
    
    index j, k, *Bdry, nB, ind[2];
    
    nB      = M->nbdry; 
    Bdry    = M->bdry;

    for ( k = 0 ; k < nB; k++)
    {
        if (Bdry[4*k+3] == 0) 
        {
            ind[0] = Bdry[4*k+0];
            ind[1] = Bdry[4*k+1];

            for (j = 0 ; j < 2 ; j++){
                u[ind[j]] = u_D;  
            }
        }
    } 


}