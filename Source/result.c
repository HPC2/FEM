#include <stdio.h>
#include <stdlib.h>
#include "hpc.h"
#include <sys/stat.h>   // stat
#include <stdbool.h>    // bool type

bool file_exists (char *filename) {
  struct stat   buffer;   
  return (stat (filename, &buffer) == 0);
}

void result_write(
    char* result_name,
    int n_rows,
    int n_cols,
    int n_refs,
    char* solver,
    int n_coords,
    int dt_build_S,
    int dt_build_rhs,
    int dt_solve
    ){
    char filename[200];
    sprintf(filename, "../Result/result_%s.csv",result_name);
    bool exists = file_exists (filename);
    FILE* f=fopen(filename, "a");
    if (!exists){
        fprintf (f,"n_rows,n_cols,n_refs,solver,n_coords,dt_build_S,dt_build_rhs,dt_solve\n");
    }
    fprintf(f,"%i,%i,%i,%s,%i,%i,%i,%i\n",n_rows,n_cols,n_refs,solver,n_coords,dt_build_S,dt_build_rhs,dt_solve);
    fclose(f);
}