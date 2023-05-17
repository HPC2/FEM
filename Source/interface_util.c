#include "hpc.h"

void interface_data_write(interface_data *interface_data, char* fname) {
    index ncoupl = interface_data->ncoupl;
    index* coupl = interface_data->coupl;
    // index* dcoupl = interface_data->dcoupl;
    // index** icoupl = interface_data->icoupl;
    // index* interface_counts = interface_data->interface_counts;
    // index** interface_ids = interface_data->interface_ids;

    char coupl_fname[64]; sprintf(coupl_fname, "%s.coupl.if", fname);

    FILE* f = fopen(coupl_fname, "w");
    for(int i = 0; i < ncoupl; i++) {
        fprintf(f, "%zu\t%zu\t%zu\t%zu\t%zu\n", 
            coupl[5*i+0],
            coupl[5*i+1],
            coupl[5*i+2],
            coupl[5*i+3],
            coupl[5*i+4]);
    }
    fclose(f);
}
