#include "hpc.h"

void print_matrix(index* A, index rows, index cols, bool rowmajor, char* name, char* extension) {
    
    index incRow = rowmajor? cols : 1;
    index incCol = rowmajor? 1 : rows;
 
    char filename[64]; sprintf(filename, "%s.%s", name, extension);

    FILE* f = fopen(filename, "w");
    for(int i=0; i<rows ; i++) {
        for(int j=0; j<cols; j++) {
            fprintf(f, "%zu\t", A[i*incRow+j*incCol]);
        }
        fprintf(f,"\n");
    }
    fclose(f);
}