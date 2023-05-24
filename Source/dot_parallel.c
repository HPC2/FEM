#include "hpc.h"
#include <mpi.h>


/*
double* dot (double* v_i, double* w_i, index n, double* ddot_res){
    ddot_res[0] = 0;
    for(int i=0; i<n; i++){
        ddot_res[0] += v_i[i] + w_i[i]; 
    }
    return ddot_res;
}
*/


double dot_parallel(double* v_i, double* w_i, index n) {
    // Input:  v_i & w_i (accumulated vectors)
    // Output: Inner product (send a scalar to all prozessors)    
    
    double glob_ddot = 0;
    double loc_ddot = 0;
    loc_ddot = ddot(v_i, 1, w_i, 1);
    MPI_ALLreduce(&loc_ddot, &glob_ddot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_ALLreduce(MPI_IN_PLACE ,&loc_ddot, 1, MPI_DOUBLE, MPI_SUM, comm);
    return glob_ddot; 

}