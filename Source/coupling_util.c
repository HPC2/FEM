#include "hpc.h"

void sort_icoupl(coupling_data* coupling){
    index   ncoupl = coupling->ncoupl;
    index*  coupl  = coupling->coupl;
    index** icoupl = coupling->icoupl;
    index*  perm[4];

    // Get memory
    index d = coupling->dcoupl[0];
    coupling->icoupl_sorted = malloc(sizeof(index*)*4);
    for(index i=0;i<4;i++){
        coupling->icoupl_sorted[i] = malloc(sizeof(index)*d);
    }
    // Copy pointer into perm
    for (int i = 0; i <4;i++){perm[i]=0;}
    for (index i = 0; i < ncoupl; i++){
        index color = coupl[5*i+4];
        perm[color] = icoupl[i];
    }
    // Fill coupl_sorted
    index* current_ptr;
    for (index i = 0; i < 4; i++){
        current_ptr = perm[i];
        for (index j = 0; j < d; j++){
            coupling->icoupl_sorted[i][j] = -1;
        }

        if (current_ptr != 0) {
            for (index j=0; j < d; j++){
                coupling->icoupl_sorted[i][j] = *(current_ptr+j);    
            }
        }
        else{
            coupling->icoupl_sorted[i]=NULL;
        }
    } 
}

void sort_coupl(coupling_data* coupling){
    index   ncoupl = coupling->ncoupl;
    index*  coupl  = coupling->coupl;
    index*  perm[4];
    
    // Get memory
    coupling->coupl_sorted = malloc(sizeof(index)*4*5);
    // Copy pointer into perm
    for (int i = 0; i <4;i++){perm[i]=0;}
    for (index i = 0; i < ncoupl; i++){
        index color = coupl[5*i+4];
        perm[color] = &coupl[5*i];
    }
    // Fill coupl_sorted 
    index* current_ptr;
    for (index i = 0; i < 4; i++){
        current_ptr = perm[i];
        if (current_ptr != 0) {
            for (index j=0; j < 5; j++){
                coupling->coupl_sorted[5*i+j] = *(current_ptr+j);    
            }
        }
        else{
            for (index j=0; j < 4; j++){
                coupling->coupl_sorted[5*i+j] = 0;    
            }
            coupling->coupl_sorted[5*i+4] = -1;
        }
    }
}

void coupling_data_print(coupling_data* coupling, int rank){
    index   n_global_cp   = coupling->n_global_cp;
    index   n_local_cp    = coupling->n_local_cp;
    index*  crossPts      = coupling->crossPts;
    index*  l2g           = coupling->l2g;
    index   ncoupl        = coupling->ncoupl;
    index*  coupl         = coupling->coupl;
    index*  coupl_sorted  = coupling->coupl_sorted;
    index*  dcoupl        = coupling->dcoupl;
    index** icoupl        = coupling->icoupl;
    index** icoupl_sorted = coupling->icoupl_sorted;

    char buf[20000];
    int offset = 0;

    // Print header and scalar data
    offset += sprintf(buf+offset,"\n");
    offset += sprintf(buf+offset,"Coupling data (rank=%i)\n", rank);
    offset += sprintf(buf+offset,"\tn_global_cp\t= %td\n", n_global_cp);
    offset += sprintf(buf+offset,"\tn_local_cp\t= %td\n", n_local_cp);
    offset += sprintf(buf+offset,"\tncoupl\t\t= %td\n", ncoupl);

    // Print cross points
    offset += sprintf(buf+offset,"\tCross points\t= [");
    for (index i=0; i<n_local_cp-1; i++){
        offset += sprintf(buf+offset,"%td, ",crossPts[i]);
    }
    offset += sprintf(buf+offset,"%td]\n", crossPts[n_local_cp-1]);

    // Print coupl
    offset += sprintf(buf+offset,"\tcoupl\t\t= [");
    for (index i=0; i<ncoupl; i++){
            offset += sprintf(
                buf+offset,
                "[%td, %td, %td, %td, %td], ",
                coupl[5*i],
                coupl[5*i+1],
                coupl[5*i+2],
                coupl[5*i+3],
                coupl[5*i+4]
            );
    }
    offset += sprintf(buf+offset,"]\n");

    // Print coupl_sorted
    offset += sprintf(buf+offset,"\tcoupl_sorted\t= [");
    for (index i=0; i<4; i++){
            offset += sprintf(
                buf+offset,
                "[%td, %td, %td, %td, %td], ",
                coupl_sorted[5*i],
                coupl_sorted[5*i+1],
                coupl_sorted[5*i+2],
                coupl_sorted[5*i+3],
                coupl_sorted[5*i+4]
            );
    }
    offset += sprintf(buf+offset,"]\n");

    // Print dcoupl
    offset += sprintf(buf+offset,"\tdcoupl\t\t= [");
    for (index i=0; i<ncoupl; i++){
            offset += sprintf(buf+offset, "%td, ", dcoupl[i]);
    }
    offset += sprintf(buf+offset,"]\n");

    // Print icoupl
    offset += sprintf(buf+offset,"\ticoupl\t\t= [");
    for (index i = 0; i < ncoupl; i++) {
        offset += sprintf(buf+offset,"[");
        for (index j = 0; j < dcoupl[i]; j++) {
            offset += sprintf(buf+offset,"%td, ", icoupl[i][j]);
        }
        offset += sprintf(buf+offset,"], ");
    }
    offset += sprintf(buf+offset,"]\n");

    // Print icoupl_sorted
    offset += sprintf(buf+offset,"\ticoupl_sorted\t= [");
    index d = dcoupl[0];
    for (index i = 0; i < 4; i++) {
        offset += sprintf(buf+offset,"[");
        if (icoupl_sorted[i]==NULL){
            offset += sprintf(buf+offset,"NULL, ");
        }
        else{
            for (index j = 0; j < d; j++) {
                offset += sprintf(buf+offset,"%td, ", icoupl_sorted[i][j]);
            }
        }
        
        offset += sprintf(buf+offset,"], ");
    }
    offset += sprintf(buf+offset,"]\n");
    printf("%s",buf);
}