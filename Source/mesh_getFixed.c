#include "hpc.h"
// Collect all nodes of typ 0 in fixedNodes 
index *mesh_getFixed(const index nCoord, const index *bdry, 
                     const index nBdry, index *nFixed)
{
    index cnt, nz, max_ind = nCoord;
    bool *flag;

    // Set a flag for each dirichlet node
    flag = (bool*) calloc(max_ind,sizeof(bool));
    cnt = 0;
    for ( unsigned int i=0 ; i < nBdry ; i++){
        if (!bdry[4*i+3]){
            cnt++;
            flag[bdry[4*i]]   = 1;
            flag[bdry[4*i+1]] = 1;
      }
    }

    // Allocate storage for fixed nodes
    index *fixed = malloc(2 * cnt * sizeof(index)); 
    if ( !(fixed) ) return(0);

    // Set the node number in fixed
    nz = 0;
    for ( unsigned int j = 0 ; j < max_ind; j++) 
    {
        if (flag[j]) fixed[nz++] = j;
    }

    *nFixed = nz;
    free(flag);
    
    return fixed;
}

