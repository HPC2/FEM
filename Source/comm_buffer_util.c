#include "hpc.h"

comm_buffers* alloc_comm_buffers(index n_global_cp, coupling_data* coupling) {
    comm_buffers* buffers = malloc(sizeof(comm_buffers));
    buffers->n_global_cp  = n_global_cp;
    buffers->cp_buffer    = calloc(n_global_cp, sizeof(double));
    buffers->if_buffer1   = calloc(coupling->dcoupl[0], sizeof(double));
    buffers->if_buffer2   = calloc(coupling->dcoupl[0], sizeof(double));
    return buffers;
}

void free_comm_buffers(comm_buffers* buffers) {
    free(buffers->cp_buffer);
    free(buffers->if_buffer1);
    free(buffers->if_buffer2);
    free(buffers);
}