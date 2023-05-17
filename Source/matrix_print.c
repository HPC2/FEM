#include "hpc.h"

void print_matrix(index* A, index rows, index cols, bool colmajor, char* name, char* extension) {
    
    index incRow = colmajor? cols : 1;
    index incCol = colmajor? 1 : rows;
 
    char filename[64]; sprintf(filename, "%s.%s", name, extension);

    FILE* f = fopen(filename, "w");
    for(int i=0; i<rows ; i++) {
        for(int j=0; j<cols; j++) {
            fprintf(f, "%zu", A[i*incRow+j*incCol]);
            if(j<cols-1) {
                fprintf(f, "\t");
            }
        }
        fprintf(f,"\n");
    }
    fclose(f);
}