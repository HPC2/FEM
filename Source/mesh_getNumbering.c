#include "hpc.h"
#include<stdio.h>

index* get_local_to_global_numbering(mesh* M, const index rows, 
                                              const index cols,
                                              index refinements) {
    index nof_nodes    = M->ncoord;
    index* elements    = M->elem;
    index nof_elements = M->nelem;

    bool* nodeslist[rows*cols];

    // allocate bool arrays for each processor
    for(int i=0; i<rows*cols; i++) {   
        nodeslist[i] = calloc(nof_nodes, sizeof(bool));
    }

    for(int i=0; i<nof_elements; i++) {
        index affiliation = elements[i*7 + 6];
        nodeslist[affiliation][elements[i*7 + 0]] = true;
        nodeslist[affiliation][elements[i*7 + 1]] = true;
        nodeslist[affiliation][elements[i*7 + 2]] = true;
    }

    index nodes_per_processor = (index) pow(pow(2,refinements)+1,2);
    printf("there are %d nodes per processor", nodes_per_processor);
    index* numbering = malloc(nodes_per_processor*(rows*cols) * sizeof(index));

    for(int i=0; i<rows*cols; i++) {
        index nodes_found = 0;
        for(int j=0; j<nof_nodes; j++) {
            if(nodeslist[i][j] == true) {
                printf("%d\t",j);
                numbering[i*nodes_per_processor + nodes_found] = j;
            }
        }
        printf("\n");
        free(nodeslist[i]);
    }

    return numbering;
}