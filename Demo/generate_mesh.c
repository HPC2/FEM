#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Pass [rows] [cols] as arguments\n");
        return -1;
    }
    int n_rows = atoi(argv[1]); //anzahl an Zeilen
    int n_cols = atoi(argv[2]); //anzahl an Spalten
    int n_proc = n_rows * n_cols;

    char pdir[] = "../Problem/";
    char coords_fname[64]; sprintf(coords_fname, "%srectangle_%dx%d.co", pdir, n_rows, n_cols);
    char elems_fname[64]; sprintf(elems_fname, "%srectangle_%dx%d.el", pdir, n_rows, n_cols);
    char bds_fname[64]; sprintf(bds_fname, "%srectangle_%dx%d.bd", pdir, n_rows, n_cols);

    // coords
    FILE *f = fopen(coords_fname, "w");
    printf("Coords:\n");
    int n_nodes_row = n_rows + 1;
    int n_nodes_col = n_cols + 1;
    int n_nodes = n_nodes_row*n_nodes_col;
    double* x_coords = (double*) malloc(sizeof(double)*n_nodes);
    double* y_coords = (double*) malloc(sizeof(double)*n_nodes);
    for (int i = 0; i < n_nodes_row; i++) {     // row_idx -> y-direction (vertical)
        for (int j = 0; j < n_nodes_col; j++) { // col_idx -> x-direction (horizontal)
            x_coords[i*n_nodes_col+j] = 1.0*j;
            y_coords[i*n_nodes_col+j] = 1.0*i;
            printf("%d:\t%d\t%d\n", i*n_nodes_col+j, j, i);
            fprintf(f, "%.1lf\t%.1lf\n", 1.0*j, 1.0*i);
        }
    }
    fclose(f);

    // elements
    f = fopen(elems_fname, "w");
    printf("Elements:\n");

    int n_elems = 2 * n_proc;
    int n_edges_h = (n_rows+1) * n_cols;
    int n_edges_v = n_rows * (n_cols + 1);

    int* elements = (int*) malloc(sizeof(int) * 7 * n_elems);
    for (int i = 0; i < n_elems; i++) {
        int alpha = i / (2 * n_cols); // column index of element
        int beta = i / 2; // affiliation
        int n0 = alpha + beta;
        int n1;
        int n2;
        int m0;
        int m1;
        int m2;
        if(i%2 == 0) {
            n1 = alpha + beta + n_cols + 2;
            n2 = alpha + beta + n_cols + 1;
            m0 = beta + n_edges_h + n_edges_v;
            m1 = beta + n_cols;
            m2 = alpha + beta + n_edges_h;
        } else {
            n1 = alpha + beta + 1;
            n2 = alpha + beta + n_cols + 2;
            m0 = beta;
            m1 = alpha + beta + n_edges_h + 1;
            m2 = beta + n_edges_h + n_edges_v;
        }
        elements[3*i+0] = n0;
        elements[3*i+1] = n1;
        elements[3*i+2] = n2;
        elements[3*i+3] = m0;
        elements[3*i+4] = m1;
        elements[3*i+5] = m2;
        elements[3*i+6] = beta;
        printf("%d:\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, n0, n1, n2, m0, m1, m2, beta);
        fprintf(f, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", n0, n1, n2, m0, m1, m2, beta);
    }
    fclose(f);

    f = fopen(bds_fname, "w");
    printf("Boundaries:\n");
    int n_edges_boundary = 2 * n_rows + 2 * n_cols;
    for (int i = 0; i < n_cols; i++) {
        int n1 = i;
        int n2 = i+1;
        int m = i;
        int type = 0;
        printf("%d\t%d\t%d\t%d\n", n1, n2, m, type);
        fprintf(f, "%d\t%d\t%d\t%d\n", n1, n2, m, type);
    }
    for (int i = 0; i < n_rows; i++) {
        int n1 = n_cols + i*n_nodes_col;
        int n2 = n_cols + (i+1)*n_nodes_col;
        int m = n_edges_h + n_cols + i*n_nodes_col;
        int type = 0;
        printf("%d\t%d\t%d\t%d\n", n1, n2, m, type);
        fprintf(f, "%d\t%d\t%d\t%d\n", n1, n2, m, type);
    }
    for (int i = 0; i < n_cols; i++) {
        int n1 = n_nodes - 1 - i;
        int n2 = n_nodes - 2 - i;
        int m = n_edges_h - 1 - i;
        int type = 0;
        printf("%d\t%d\t%d\t%d\n", n1, n2, m, type);
        fprintf(f, "%d\t%d\t%d\t%d\n", n1, n2, m, type);
    }
    for (int i = 0; i < n_rows; i++) {
        int n1 = n_nodes - 1 - (n_cols + i*n_nodes_col);
        int n2 = n_nodes - 1 - (n_cols + (i+1)*n_nodes_col);
        int m = n_edges_h - 1 + n_edges_v - (n_cols + i*n_nodes_col);
        int type = 0;
        printf("%d\t%d\t%d\t%d\n", n1, n2, m, type);
        fprintf(f, "%d\t%d\t%d\t%d\n", n1, n2, m, type);
    }
    fclose(f);


    free(x_coords);
    free(y_coords);
    free(elements);
}