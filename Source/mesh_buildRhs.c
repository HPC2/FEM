#include "hpc.h"

// Boundary integrals (g, \varphi)_L2(\Gamma_N)
void mesh_Neumann(double p1[2], double p2[2], 
                  index typ, 
                  double (*fc)(double *, index), 
                  double m[2])  
{
    int i;
    double mid[2], edge_len = 0.0, fac;

    // Get the length of the edge
    for (i = 0 ; i < 2 ; i++ ){
        edge_len += (p2[i] - p1[i]) * (p2[i] - p1[i]) ;
        mid[i]    = (p1[i] + p2[i])/2.0;
    }
    edge_len = sqrt(edge_len);

    // Eval midpoint rule
    fac     = (1./2.)*fc(mid,typ)*edge_len;

    // Assign values
    m[0]    = fac; 
    m[1]    = fac;
}

// Calculate the values for the integral (f,\varphi)_T
// ordering w.r.t. [ p1, p2, p3] 
void mesh_vol_elem(double p1[2], double p2[2], double p3[2], 
              index typ, double (*fc)(double *, index), double m[3])  
{
    int i;
    double mid[2], d[2][2], fac;

    // Get vectors and midpoint
    for (i = 0 ; i < 2 ; i++ ){
        d[0][i] = p1[i]-p3[i];
        d[1][i] = p2[i]-p1[i];
        mid[i]  = (p1[i] + p2[i] + p3[i])/3.0;
    }

    // Get value with midpoint rule
    fac = fc(mid,typ)/6.0*(d[0][0]*d[1][1]-d[1][0]*d[0][1]);

    // Every node has the same value
    for ( i = 0 ; i < 3 ; i++ ){
        m[i] = fac; 
    }
}


// Compute right hand side using midpoint rule
void mesh_buildRhs(const mesh *M, double *b, double (*fV)(double *, index), 
                   double (*fN)(double *, index))
{  
    index j, k, nC, nT, nB, *Elem, *Bdry, ind[3] ;
    double *Coord, bx[3] ;

    // For convenience
    nC      = M->ncoord; 
    nT      = M->nelem; 
    nB      = M->nbdry; 
    Coord   = M->coord; 
    Elem    = M->elem; 
    Bdry    = M->bdry; 
    
    // Calculate the values for the integral (f,\varphi)_\Omega
    for ( k = 0 ; k < nT; k++)
    {
        // Get the element nodes
        for (j = 0 ; j < 3 ; j++){
            ind[j] = Elem[7*k+j];
        }

        // Calculate the integral over the element
        mesh_vol_elem(Coord+2*ind[0],
                      Coord+2*ind[1],
                      Coord+2*ind[2],
                      Elem[7*k+6],
                      fV,
                      bx);

        // Add up the values to global rhs vector
        for (j = 0 ; j < 3 ; j++){
            b[ind[j]] += bx[j] ;  
        }
    }                    
    
    // Boundary integrals (g, \varphi)_L2(\Gamma_N)
    for ( k = 0 ; k < nB; k++)
    {
        if (Bdry[4*k+3]) // if neumann
        {
            // Get edge nodes
            ind[0] = Bdry[4*k+0];
            ind[1] = Bdry[4*k+1];

            // Calculate integral over edge
            mesh_Neumann(Coord+2*ind[0],Coord+2*ind[1],Bdry[4*k+3],fN,bx);

            // Add up the values to global rhs vector
            for (j = 0 ; j < 2 ; j++){
                b[ind[j]] += bx[j];  
            }
        }
    }                    
}



