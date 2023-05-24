#include "hpc.h"

void write_l2g(index* l2g_numbering, index n_rows, index n_cols, index refinements) {

    char fname[200];
    sprintf(fname,"../Problem/l2g_%tdx%td_%tdref",n_rows,n_cols,refinements);
    print_imatrix(l2g_numbering, n_rows*n_cols, (index) pow(pow(2,refinements)+1,2) , true, fname, "l2g");
}