#include "hpc.h"
#include <mpi.h>

index* mpi_boundaries(index n_rows, index n_cols, index* global_boundaries) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_p; MPI_Comm_size(MPI_COMM_WORLD, &nof_p);
    index* boundaries_;
    index* boundaries = malloc(sizeof(index)*4);
    if (rank == 0) {
        index n_proc = n_rows * n_cols;
        boundaries_ = malloc(sizeof(index)*n_proc*4); // boundaries for each proc
        for (index i = 0; i < n_proc; i++) {
            boundaries_[i*4 + 0] = i < n_cols ? global_boundaries[0] : -1;
            boundaries_[i*4 + 1] = (i+1) % n_cols == 0 ? global_boundaries[1] : -1;
            boundaries_[i*4 + 2] = i >= n_proc - n_cols ? global_boundaries[2] : -1;
            boundaries_[i*4 + 3] = i % n_cols == 0 ? global_boundaries[3] : -1;
        } 
    }
    MPI_Datatype bd_type;
    MPI_Type_contiguous(4, MPI_AINT, &bd_type);
    MPI_Type_commit(&bd_type);
    MPI_Scatter(boundaries_, 1, bd_type, boundaries, 1, bd_type, 0, MPI_COMM_WORLD);
    MPI_Type_free(&bd_type);
    return boundaries;

}