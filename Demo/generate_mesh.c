#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Pass [rows] [cols] as arguments\n");
        return -1;
    }
    int n_p_rows = atoi(argv[1]); //anzahl an Zeilen
    int n_p_cols = atoi(argv[2]); //anzahl an Spalten
    int n_p = n_p_rows * n_p_cols;

    // coords
    FILE *f = fopen("custom_problem.co", "w");
    printf("Coords:\n");
    int n_nodes_row = n_p_rows + 1;
    int n_nodes_col = n_p_cols + 1;
    int n_nodes = n_nodes_row*n_nodes_col;
    double* x_coords = malloc(sizeof(double)*n_nodes);
    double* y_coords = malloc(sizeof(double)*n_nodes);
    for (int i = 0; i < n_nodes_row; i++) {
        for (int j = 0; j < n_nodes_col; j++) {
            x_coords[i*n_nodes_col+j] = 1.0*i;
            y_coords[i*n_nodes_col+j] = 1.0*j;
            printf("%d:\t%d\t%d\n", i*n_nodes_col+j, i, j);
            fprintf(f, "%.1lf\t%.1lf\n", 1.0*i, 1.0*j);
        }
    }
    fclose(f);

    // elements
    f = fopen("custom_problem.el", "w");
    printf("Elements:\n");
    int n_elem = 2*n_p;
    int* elements = malloc(sizeof(int) * 3 * n_elem);
    for (int e = 0; e < n_elem; e++) {
        int n0;
        int n1;
        int n2;
        if(e%2 == 0) {
            n0 = e/2 + e/(2*n_p_cols);
            n1 = n0 + n_nodes_col + 1;
            n2 = n1 - 1;
        } else {
            n0 = (e-1)/2 + (e-1)/(2*n_p_cols);
            n1 = n0 + 1;
            n2 = n0 + n_nodes_col + 1;
        }
        elements[3*e+0] = n0;
        elements[3*e+1] = n1;
        elements[3*e+2] = n2;
        printf("%d:\t%d\t%d\t%d\n", e, n0, n1, n2);
        fprintf(f, "%d\t%d\t%d\n", n0, n1, n2);
    }

    int n_diag_edges = n_p;
    int n_horz_edges = (n_p_rows+1) * n_p_cols;
    int n_vert_edges = n_p_rows * (n_p_cols + 1);
    
    // TODO: Edges

    fclose(f);




    free(x_coords);
    free(y_coords);
}