#include "hpc.h"
#include<stdio.h>
#include<assert.h>

index* get_local_to_global_numbering(mesh* M, const index rows, 
                                              const index cols,
                                              index refinements) {
/*
returns an index matrix where the rows are the permutation vector for each processor
*/

    index nof_nodes    = M->ncoord;
    printf("there are %zu nodes in total\n", nof_nodes);
    index* elements    = M->elem;
    index nof_elements = M->nelem;

    bool* nodeslist[rows*cols];

    // allocate bool arrays for each processor
    for(int i=0; i<rows*cols; i++) {   
        nodeslist[i] = (bool*) calloc(nof_nodes, sizeof(bool)); //zero-allocated
    }

    // set node numbers of each processor/affiliation to true
    for(int i=0; i<nof_elements; i++) {
        index affiliation = elements[i*7 + 6];
        //printf("processor %zu: found nodes %zu, %zu and %zu\n", affiliation,elements[i*7 + 0],elements[i*7 + 1],elements[i*7 + 2]);
        nodeslist[affiliation][elements[i*7 + 0]] = true;
        nodeslist[affiliation][elements[i*7 + 1]] = true;
        nodeslist[affiliation][elements[i*7 + 2]] = true;
        //debugging
        /*
        for(int i=0; i<rows*cols; i++) {
        index nodes_found = 0;
            for(int j=0; j<nof_nodes; j++) {
                printf("%d", nodeslist[i][j]);
            }
            printf("\n");
        }
        */
    }

    index nodes_per_processor = (index) pow(pow(2,refinements)+1,2); // (2^f + 1)^2
    printf("there are %zu nodes per processor:\n", nodes_per_processor);
    index* numbering = malloc(nodes_per_processor*(rows*cols) * sizeof(index));

    for(int i=0; i<rows*cols; i++) {
        index nodes_found = 0;
        for(int j=0; j<nof_nodes; j++) {
            //printf("%d", nodeslist[i][j]);
            if(nodeslist[i][j] == true) {
                //printf("%d ",j);
                numbering[i*nodes_per_processor + nodes_found++] = j;
            }
        }
        assert(nodes_found == nodes_per_processor);
        //printf("\n");
    }

    for(int i=0; i<rows*cols; i++) {
        free(nodeslist[i]);
    }

    return numbering;
}

index* get_global_to_local_numbering(index* local_numbering, 
                                     index  nof_local_nodes, // nodes per processor = (2^f+1)^2
                                     index  nof_global_nodes) {
/*
Inputs: -local-to-global numbering of ONE processor
        -number of local nodes
        -number of global nodes

returns the permutation vector for the global to local numbering
non-existent nodes are marked with -1
*/                  
    index* numbering = (index*) malloc(nof_global_nodes*sizeof(index));
    for(int j=0; j<nof_global_nodes; j++) {
        numbering[j] = -1;
    }

    for(int i=0; i<nof_local_nodes; i++) {
        numbering[local_numbering[i]] = i;
    }

    return numbering;
}