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
    int n_processors,
    int n_rows,
    int n_cols,
    int n_refs,
    char* solver,
    int n_iter_max,
    int n_iter_min,
    int n_iter_sum,
    int n_coords,
    index* boundaries,
    size_t dt_init,
    size_t dt_build_S,
    size_t dt_build_rhs,
    size_t dt_solve,
    size_t t_cp_comm,
    size_t t_if_comm
    ){
    char filename[200];
    sprintf(filename, "../Result/result_%s.csv",result_name);
    bool exists = file_exists (filename);
    FILE* f=fopen(filename, "a");
    if (!exists){
        fprintf (f,"n_rows,");
        fprintf (f,"n_cols,");
        fprintf (f,"n_refs,");
        fprintf (f,"solver,");
        fprintf (f,"n_iter_max,");
        fprintf (f,"n_iter_min,");
        fprintf (f,"n_iter_sum,");
        fprintf (f,"n_processors,");
        fprintf (f,"n_coords,");
        fprintf (f,"dt_init,");
        fprintf (f,"dt_build_S,");
        fprintf (f,"dt_build_rhs,");
        fprintf (f,"dt_solve,");
        fprintf (f,"t_cp_comm,");
        fprintf (f,"t_if_comm,");
        fprintf (f,"glob_b0,");
        fprintf (f,"glob_b1,");
        fprintf (f,"glob_b2,");
        fprintf (f,"glob_b3\n");
        
    }
    fprintf(f,"%i,%i,%i,%s,%i,%i,%i,%i,%i,%zu,%zu,%zu,%zu,%zu,%zu,%td,%td,%td,%td\n",
      n_rows,
      n_cols,
      n_refs,
      solver,
      n_iter_max,
      n_iter_min,
      n_iter_sum,
      n_processors,
      n_coords,
      dt_init,
      dt_build_S,
      dt_build_rhs,
      dt_solve,
      t_cp_comm,
      t_if_comm,
      boundaries[0],
      boundaries[1],
      boundaries[2],
      boundaries[3]
      );
    fclose(f);
}