#include "hpc.h"
#include<stdio.h>

index* get_local_to_global_numbering(mesh* M, const index rows, 
                                              const index cols) {
    index nof_nodes = M->ncoord;
    index* elements = M->elem;

    bool* nodeslist[rows*cols];

    for(int i=0; i<rows*cols; i++) {   
        nodeslist[i] = calloc(nof_nodes, sizeof(bool));
    }

    for(int i=0; i<rows*cols; i++) {
        index affiliation = elements[i*7 + 6];
        nodeslist[affiliation][elements[i*7 + 0]] = true;
        nodeslist[affiliation][elements[i*7 + 1]] = true;
        nodeslist[affiliation][elements[i*7 + 2]] = true;
    }

    for(int i=0; i<rows*cols; i++) {
        for(int j=0; nof_nodes; j++) {
            if(nodeslist[i][j] == true) {
                printf("%d ",j);
            }
        }
        printf("\n");
        free(nodeslist[i]);
    }
}