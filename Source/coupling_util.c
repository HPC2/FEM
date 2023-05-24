#include "hpc.h"

void coupling_data_print(coupling_data* coupling, int rank){
    index   nCrossPts = coupling->n_local_cp;
    index*  crossPts  = coupling->crossPts;
    index*  l2g       = coupling->l2g;
    index   ncoupl    = coupling->ncoupl;
    index*  coupl     = coupling->coupl;
    index*  dcoupl    = coupling->dcoupl;
    index** icoupl    = coupling->icoupl;

    char buf[20000];
    int offset = 0;

    // Print header and scalar data
    offset += sprintf(buf+offset,"\n");
    offset += sprintf(buf+offset,"Coupling data (rank=%i)\n", rank);
    offset += sprintf(buf+offset,"\tnCrossPts\t= %d\n", nCrossPts);
    offset += sprintf(buf+offset,"\tncoupl\t\t= %d\n", ncoupl);

    // Print cross points
    offset += sprintf(buf+offset,"\tCross points\t= [");
    for (index i=0; i<nCrossPts-1; i++){
        offset += sprintf(buf+offset,"%td, ",crossPts[i]);
    }
    offset += sprintf(buf+offset,"%d]\n", crossPts[nCrossPts-1]);

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
    
    printf("%s",buf);
}