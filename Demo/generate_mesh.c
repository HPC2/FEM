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
    
    // refine the mesh
    int refinements = 1;
    for(int i=0; i<refinements; i++) {
        Mesh = mesh_refine(Mesh);
        // update edge/node info??
        mesh_getEdge2no( Mesh->nelem, 
                     Mesh->elem, 
                    &Mesh->nedges,
                    &Mesh->edge2no);
    }

    mesh_print(Mesh, 0);
    char fname[64]; sprintf(fname, "%srectangle_%dx%d_refined%dtimes", pdir, n_rows, n_cols, refinements);
    mesh_write(Mesh, fname);

    mesh_free(Mesh);
}