#include <stdio.h>
#include <stdlib.h>
#include "hpc.h"

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Pass [rows] [cols] as arguments\n");
        return -1;
    }

    int n_rows = atoi(argv[1]); //anzahl an Zeilen
    int n_cols = atoi(argv[2]); //anzahl an Spalten

    char pdir[] = "../Problem/";

    bool boundaries[4] = {0, 0, 0, 0};
    mesh* Mesh = mesh_create_rect(n_rows, n_cols, boundaries);
    if(!Mesh) {
        printf("OOM\n");
        return 1;
    }

    mesh_getEdge2no( Mesh->nelem, 
                     Mesh->elem, 
                    &Mesh->nedges,
                    &Mesh->edge2no);
    // mesh_print(Mesh, 0);
    char fname[64]; sprintf(fname, "%srectangle_%dx%d_refined", pdir, n_rows, n_cols);


    Mesh = mesh_refine(Mesh);

    mesh_getEdge2no( Mesh->nelem, 
                     Mesh->elem, 
                    &Mesh->nedges,
                    &Mesh->edge2no);

    interface_data* interfaces = rect_interface_data(Mesh, n_rows, n_cols, 1);
    for (int i = 0; i < interfaces->ncoupl; i++) {
        printf("%zu\n", interfaces->interf2edge[i]);
        for (int j = 0; j < 5; j++) {
            printf("%zu\t", interfaces->coupl[5*i+j]);
        }
        printf("\n");
        for (int j = 0; j < interfaces->dcoupl[i]; j++) {
            printf("%zu\t", interfaces->icoupl[i][j]);
        }
        printf("\n\n");
    }
    // mesh_write(Mesh, fname);
    mesh_free(Mesh);
}