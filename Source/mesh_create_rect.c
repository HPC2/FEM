#include "hpc.h"
#include <stdio.h>

mesh *mesh_create_rect(index n_rows, index n_cols, bool *boundaries) {
    index n_nodes_row = n_rows + 1;
    index n_nodes_col = n_cols + 1;
    index n_nodes = n_nodes_row*n_nodes_col;
    index n_elems = 2 * n_rows * n_cols;
    index n_edges_boundary = 2 * n_rows + 2 * n_cols;

    mesh* Mesh = mesh_alloc(n_nodes, n_elems, n_edges_boundary);

    // coords
    double* coords = Mesh->coord;
    for (index i = 0; i < n_nodes_row; i++) {     // row_idx -> y-direction (vertical)
        for (index j = 0; j < n_nodes_col; j++) { // col_idx -> x-direction (horizontal)
            index node_id = i*n_nodes_col+j;
            coords[2*node_id+0] = 1.0*j;
            coords[2*node_id+1] = 1.0*i;
        }
    }

    //elements
    index n_edges_h = (n_rows+1) * n_cols;      // num horizontal edges
    index n_edges_v = n_rows * (n_cols + 1);    // num vertical edges
    index* elements = Mesh->elem;
    for (index i = 0; i < n_elems; i++) {
        index alpha = i / (2 * n_cols); // column index of element
        index beta = i / 2; // affiliation
        index n0 = alpha + beta;
        index n1;
        index n2;
        index m0;
        index m1;
        index m2;
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
        elements[7*i+0] = n0;
        elements[7*i+1] = n1;
        elements[7*i+2] = n2;
        elements[7*i+3] = m0;
        elements[7*i+4] = m1;
        elements[7*i+5] = m2;
        elements[7*i+6] = beta;
    }

    // boundary edges
    index* bdry = Mesh->bdry;
    index offset = 0;
    for (index i = 0; i < n_cols; i++) {
        bdry[offset++] = i;
        bdry[offset++] = i+1;
        bdry[offset++] = i;
        bdry[offset++] = boundaries[0];
    }
    for (index i = 0; i < n_rows; i++) {
        bdry[offset++] = n_cols + i*n_nodes_col;
        bdry[offset++] = n_cols + (i+1)*n_nodes_col;
        bdry[offset++] = n_edges_h + n_cols + i*n_nodes_col;
        bdry[offset++] = boundaries[1];
    }
    for (index i = 0; i < n_cols; i++) {
        bdry[offset++] = n_nodes - 1 - i;
        bdry[offset++] = n_nodes - 2 - i;
        bdry[offset++] = n_edges_h - 1 - i;
        bdry[offset++] = boundaries[2];
    }
    for (index i = 0; i < n_rows; i++) {
        bdry[offset++] = n_nodes - 1 - (n_cols + i*n_nodes_col);
        bdry[offset++] = n_nodes - 1 - (n_cols + (i+1)*n_nodes_col);
        bdry[offset++] = n_edges_h - 1 + n_edges_v - (n_cols + i*n_nodes_col);
        bdry[offset++] = boundaries[3];
    }

    return Mesh;
}

interface_data* rect_interface_data(mesh *Mesh, index n_rows, index n_cols, index n_refinments) {
    index n_interfaces_h = (n_rows-1)*n_cols;
    index n_interfaces_v = (n_cols-1)*n_rows;
    interface_data* interfaces = malloc(sizeof(interface_data));
    // crosspoints
    interfaces->nCrossPts = Mesh->ncoord;
    interfaces->crossPts = malloc(sizeof(index) * interfaces->nCrossPts);
    for (int i = 0; i < interfaces->nCrossPts; i++) {
        interfaces->crossPts[i] = i;
    }
    // interfaces
    index n_interfaces = n_interfaces_h + n_interfaces_v;
    interfaces->ncoupl = n_interfaces;
    interfaces->coupl = malloc(sizeof(index) * 5 * n_interfaces);
    interfaces->interf2edge = malloc(sizeof(index) * n_interfaces);
    index* coupl = interfaces->coupl;
    index offset = 0;
    for (index i = 0; i < n_interfaces_h; i++) {
        index beta = i / n_cols;
        coupl[offset++] = beta + n_cols + i + 1;
        coupl[offset++] = beta + n_cols + i + 2;
        coupl[offset++] = n_cols + i;
        coupl[offset++] = i;
        coupl[offset++] = 0; // TODO: Color
        interfaces->interf2edge[i] = n_cols + i;
    }
    for (index i = 0; i < n_interfaces_v; i++) {
        index alpha = i % n_rows;
        index beta = i / n_rows;
        index gamma = i % (n_cols-1);
        index delta = i / (n_cols-1);
        coupl[offset++] = (n_cols+1)*alpha + beta +1;
        coupl[offset++] = (n_cols+1)*alpha + beta +1 + n_cols + 1;
        coupl[offset++] = n_cols*alpha + beta;
        coupl[offset++] = n_cols*alpha + beta + 1;
        coupl[offset++] = 0; // TODO: Color
        interfaces->interf2edge[n_interfaces_h+i] = gamma + 
            delta * (n_cols + 1) + 
            n_cols*(n_rows+1) + 1;
    }


    // interface nodes
    index* edge2no = Mesh->edge2no;
    size_t nnpi = 1 << n_refinments - 1; // num (inner) nodes per interface
    interfaces->icoupl = malloc(sizeof(index*) * n_interfaces);
    interfaces->dcoupl = malloc(sizeof(index) * n_interfaces);
    // printf("nnpi:%d\n", nnpi);
    // printf("nedges:%d\n", Mesh->nedges);
    for (index i = 0; i < n_interfaces; i++) {
        interfaces->dcoupl[i] = nnpi;
        index original_edge = interfaces->interf2edge[i];
        interfaces->icoupl[i] = malloc(sizeof(index) * nnpi);
        index* icoupl = interfaces->icoupl[i];
        printf("i:%d\n", i);
        for (index j = 0; j < nnpi; j++) {
            printf("%zu, %zu\n", edge2no[2*((nnpi+1)*original_edge + j) + 0], edge2no[2*((nnpi+1)*original_edge + j) + 1]);
            icoupl[j] = edge2no[2*((nnpi+1)*original_edge + j) + 1];
        }
        
    }

    return interfaces;
}