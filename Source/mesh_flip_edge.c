#include "hpc.h"
#include <stdio.h>

void mesh_flip_edge(mesh* m) {
    index* elem = m->elem;
    index nelem = m->nelem;
    index lr_elem_id = -1;
    index ul_elem_id = -1;
    index lower_right_elem[7];
    index upper_left_elem[7];
    for (index i = 0; i < nelem; i++) {
        if (elem[i*7] == 1) {
            lr_elem_id = i;
            for (int j = 0; j < 7; j++) {
                lower_right_elem[j] = elem[i*7+j];
            }
        } else if (elem[i*7] == 2) {
            ul_elem_id = i;
            for (int j = 0; j < 7; j++) {
                upper_left_elem[j] = elem[i*7+j];
            }
        }
    }
    // printf("lower right elem: %td\n", lr_elem_id);
    // for (int i = 0; i < 6; i++) {
    //     printf("%td\t", lower_right_elem[i]);
    // }
    // printf("\n");
    // printf("upper left elem: %td\n", ul_elem_id);
    // for (int i = 0; i < 6; i++) {
    //     printf("%td\t", upper_left_elem[i]);
    // }
    // printf("\n");

    index l_edge_replace = lower_right_elem[4];
    //printf("edge to replace: %td\n", l_edge_replace);
    index u_edge_replace = upper_left_elem[4];
    //printf("edge to replace: %td\n", u_edge_replace);

    index l_other_elem_id = -1;
    index u_other_elem_id = -1;
    index l_other_elem[7];
    index u_other_elem[7];
    for (index i = 0; i < nelem; i++) {
        for (index j = 3; j < 6; j++) {
            if (elem[i*7+j] == l_edge_replace) {
                l_other_elem_id = i;
                for (int k = 0; k < 7; k++) {
                    l_other_elem[k] = elem[i*7+k];
                }
            } else if (elem[i*7+j] == u_edge_replace) {
                u_other_elem_id = i;
                for (int k = 0; k < 7; k++) {
                    u_other_elem[k] = elem[i*7+k];
                }
            }
        }
    }

    // printf("lower right other elem: %td\n", l_other_elem_id);
    // for (int i = 0; i < 6; i++) {
    //     printf("%td\t", l_other_elem[i]);
    // }
    // printf("\n");
    // printf("upper left other elem: %td\n", u_other_elem_id);
    // for (int i = 0; i < 6; i++) {
    //     printf("%td\t", u_other_elem[i]);
    // }
    // printf("\n");

    index l_other_node = -1;
    index l_other_node_pos = -1;
    for (int i = 0; i < 3; i++) {
        if (l_other_elem[i] != lower_right_elem[1] &&
            l_other_elem[i] != lower_right_elem[2]) {
            l_other_node = l_other_elem[i];
            l_other_node_pos = i;
        }
    }

    // printf("lower right other node: %td\n", l_other_node);
    index nodes[4] = {
        lower_right_elem[0], 
        lower_right_elem[1], 
        l_other_node, 
        lower_right_elem[2]
    };

    index edges[4] =  {
        lower_right_elem[3],
        l_other_elem[3 + (l_other_node_pos+2)%3],
        l_other_elem[3 + l_other_node_pos],
        lower_right_elem[5]
    };

    //for (int i = 0; i < 4; i++){
    //    printf("%td\t", nodes[i]);
    //}
    //printf("\n");

    //for (int i = 0; i < 4; i++){
    //    printf("%td\t", edges[i]);
    //}
    //printf("\n");

    elem[7*lr_elem_id+0] = nodes[0];
    elem[7*lr_elem_id+1] = nodes[1];
    elem[7*lr_elem_id+2] = nodes[2];
    elem[7*lr_elem_id+3] = edges[0];
    elem[7*lr_elem_id+4] = edges[1];
    elem[7*lr_elem_id+5] = l_edge_replace;

    elem[7*l_other_elem_id+0] = nodes[2];
    elem[7*l_other_elem_id+1] = nodes[3];
    elem[7*l_other_elem_id+2] = nodes[0];
    elem[7*l_other_elem_id+3] = edges[2];
    elem[7*l_other_elem_id+4] = edges[3];
    elem[7*l_other_elem_id+5] = l_edge_replace;


    index u_other_node = -1;
    index u_other_node_pos = -1;
    for (int i = 0; i < 3; i++) {
        if (u_other_elem[i] != upper_left_elem[1] &&
            u_other_elem[i] != upper_left_elem[2]) {
            u_other_node = u_other_elem[i];
            u_other_node_pos = i;
        }
    }
    
    // printf("upper left other node: %td\n", u_other_node);
    index u_nodes[4] = {
        upper_left_elem[0], 
        upper_left_elem[1], 
        u_other_node, 
        upper_left_elem[2]
    };

    index u_edges[4] =  {
        upper_left_elem[3],
        u_other_elem[3 + (u_other_node_pos+2)%3],
        u_other_elem[3 + u_other_node_pos],
        upper_left_elem[5]
    };

    //for (int i = 0; i < 4; i++){
    //    printf("%td\t", u_nodes[i]);
    //}
    //printf("\n");

    //for (int i = 0; i < 4; i++){
    //    printf("%td\t", u_edges[i]);
    //}
    //printf("\n");

    elem[7*ul_elem_id+0] = u_nodes[0];
    elem[7*ul_elem_id+1] = u_nodes[1];
    elem[7*ul_elem_id+2] = u_nodes[2];
    elem[7*ul_elem_id+3] = u_edges[0];
    elem[7*ul_elem_id+4] = u_edges[1];
    elem[7*ul_elem_id+5] = u_edge_replace;

    elem[7*u_other_elem_id+0] = u_nodes[2];
    elem[7*u_other_elem_id+1] = u_nodes[3];
    elem[7*u_other_elem_id+2] = u_nodes[0];
    elem[7*u_other_elem_id+3] = u_edges[2];
    elem[7*u_other_elem_id+4] = u_edges[3];
    elem[7*u_other_elem_id+5] = u_edge_replace;

    mesh_getEdge2no( m->nelem, 
                     m->elem, 
                    &m->nedges,
                    &m->edge2no);
}