#include "hpc.h"

// Refine uniformly a regular mesh
mesh *mesh_refine(const mesh *In)
{
    // Declare new and aux variables (n = new)
    index i, i0, i1, j, k, p, nCoord, nElem, nEdges, nBdry, *E, *nE, 
            *B, *nB, *edge2no;
    // Aux indeces
    static index isucc[3] = {1,2,0}, iprae[3] = {2,0,1};
    // Old and new coordinates
    double *C, *nC;
    // Refined mesh
    mesh *Out;

    // Get the data structures from initial mesh
    C = In->coord; 
    E = In->elem; 
    B = In->bdry;
    nCoord = In->ncoord; 
    nElem  = In->nelem; 
    nBdry  = In->nbdry;    
    nEdges = In->nedges;

    // Allocate storage for refined mesh */
    Out = mesh_alloc(nCoord+nEdges, 4*nElem, 2*nBdry);

    // Get the old edge to node information
    edge2no = In->edge2no;

    // Declare convinient variables
    nC = Out->coord; 
    nE = Out->elem; 
    nB = Out->bdry;

    // Copy old coordinates
    for ( i = 0 ; i < 2*nCoord ; i++) 
    {
        nC[i] = C[i];
    }

    // Compute new coordinates
    for ( i = 0 ; i < nEdges ; i++)
    {
        i0 = edge2no[2*i]  ; 
        i1 = edge2no[2*i+1];  
        nC[2*nCoord+2*i  ] = 0.5*(C[2*i0  ] + C[2*i1  ]);
        nC[2*nCoord+2*i+1] = 0.5*(C[2*i0+1] + C[2*i1+1]);
    }   

    // Create new elements
    for ( i = 0 ; i < nElem ; i++)
    {
        for ( k = 0 ; k < 3 ; k++)
        {
            // Element defining nodes
            nE[28*i+7*k+0] = E[7*i+k];
            nE[28*i+7*k+1] = nCoord + E[7*i+3+k];
            nE[28*i+7*k+2] = nCoord + E[7*i+3+iprae[k]];

            // Edge number
            nE[28*i+7*k+3] = 2*E[7*i+3+k]        + (E[7*i+k]>E[7*i+isucc[k]]);
            nE[28*i+7*k+4] = 2*nEdges+3*i+k;
            nE[28*i+7*k+5] = 2*E[7*i+3+iprae[k]] + (E[7*i+k]>E[7*i+iprae[k]]);

            // Affiliation of element
            nE[28*i+7*k+6] = E[7*i+6];         
        }
        nE[28*i+21] = nCoord + E[7*i+3];
        nE[28*i+22] = nCoord + E[7*i+4];
        nE[28*i+23] = nCoord + E[7*i+5];

        nE[28*i+24] = 2*nEdges+3*i+1;
        nE[28*i+25] = 2*nEdges+3*i+2;
        nE[28*i+26] = 2*nEdges+3*i+0;

        nE[28*i+27] = E[7*i+6];         
    }  

    // Create refined boundary
    for ( i = 0 ; i < nBdry ; i++ ){
        nB[8*i  ] = B[4*i];  
        nB[8*i+1] = nCoord+B[4*i+2];
        nB[8*i+4] = nCoord+B[4*i+2];
        nB[8*i+5] = B[4*i+1];
            
        nB[8*i+2] = 2*B[4*i+2] + (B[4*i  ] > B[4*i+1]) ;
        nB[8*i+6] = 2*B[4*i+2] + (B[4*i+1] > B[4*i  ]) ;
            
        nB[8*i+3] = B[4*i+3] ;
        nB[8*i+7] = B[4*i+3]  ;
    }

    return (Out) ;
}



