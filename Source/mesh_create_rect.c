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
    index n_proc = n_rows*n_cols;
    interface_data* interfaces = malloc(sizeof(interface_data));

    // interfaces
    index n_interfaces = n_interfaces_h + n_interfaces_v;
    interfaces->ncoupl = n_interfaces;
    interfaces->coupl = malloc(sizeof(index) * 5 * n_interfaces);
    index* interf2edge = malloc(sizeof(index) * n_interfaces);
    index* coupl = interfaces->coupl;
    index offset = 0;
    index source_offset        = n_cols + 1;
    index target_offset        = n_cols + 2;
    index left_process_offset  = n_cols;
    index right_process_offset = 0;
    index midpoint_offset      = n_cols;
    for (index i = 0; i < n_interfaces_h; i++) {
        index alpha = i % n_cols;
        index beta = i / n_cols;
        coupl[offset++] = beta + i + source_offset;
        coupl[offset++] = beta + i + target_offset;
        coupl[offset++] = i + left_process_offset;
        coupl[offset++] = i + right_process_offset;
        coupl[offset++] = (alpha + (beta % 2)) % 2;
        interf2edge[i] = i + midpoint_offset;
    }
    source_offset        = 1;
    target_offset        = n_cols + 2;
    left_process_offset  = 0;
    right_process_offset = 1;
    midpoint_offset      = n_cols*(n_rows+1) + 1;
    for (index i = 0; i < n_interfaces_v; i++) {
        index alpha = i % (n_cols-1);
        index beta = i/(n_cols-1);
        coupl[offset++] = 2 * beta + i + source_offset;
        coupl[offset++] = 2 * beta + i + target_offset;
        coupl[offset++] = beta + i + left_process_offset;
        coupl[offset++] = beta + i + right_process_offset;
        coupl[offset++] = (alpha + (beta % 2) + 1) % 2 + 2;
        interf2edge[n_interfaces_h+i] = 
            2 * beta + i + midpoint_offset;
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
        index original_edge = interf2edge[i];
        interfaces->icoupl[i] = malloc(sizeof(index) * nnpi);
        index* icoupl = interfaces->icoupl[i];
        // printf("i:%d\n", i);
        for (index j = 0; j < nnpi; j++) {
            // printf("%zu, %zu\n", edge2no[2*((nnpi+1)*original_edge + j) + 0], edge2no[2*((nnpi+1)*original_edge + j) + 1]);
            icoupl[j] = edge2no[2*((nnpi+1)*original_edge + j) + 1];
        }
        
    }

    index* interface_counts = calloc(n_proc, sizeof(index));
    interfaces->interface_ids = malloc(sizeof(index*)*n_proc);
    index** interface_ids = interfaces->interface_ids;
    for (int i = 0; i < n_proc; i++) {
        interface_ids[i] = malloc(sizeof(index)*4);
    }
    for (int i = 0; i < interfaces->ncoupl; i++) {
        index left = interfaces->coupl[i*5 + 2];
        index right = interfaces->coupl[i*5 + 3];
        interface_ids[left][interface_counts[left]] = i;
        interface_ids[right][interface_counts[right]] = i;
        interface_counts[left]++;
        interface_counts[right]++;
    }
    interfaces->interface_counts = interface_counts;

    return interfaces;
}