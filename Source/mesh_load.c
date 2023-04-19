#include "hpc.h"

// Load a (index) matrix from a file (usually elements or bdry data)
index *mesh_load_index (char *fname, index cols, index *rows)
{
    FILE *file;
    index cnt, j, a, *data;
            
    // try open file
    file = fopen(fname,"r");
    // exception handling
    if (file == NULL) return (NULL) ;

    // count entries
    cnt = 0;
    while (fscanf(file,"%zu",&a) != EOF) cnt++;

    // check if matrix is square
    rows[0] = cnt/cols;
    if (  cnt - rows[0] * cols ){
      fclose(file);
      printf("\nmesh_load_index() Error!!! cnt = %g, rows %g\n\n", 
              (double) cnt, (double) rows[0]);
      return (NULL);
    }
    
    // allocate memory
    data = malloc(cnt * sizeof(index));
    if (!data) return (NULL) ;
    
    // transmit data to data
    fseek(file,0L,SEEK_SET);
    for (j=0; j<cnt; j++) {
      fscanf(file,"%zu",&(data[j]));
    }

    // close file
    fclose(file);
    return (data);
}

// Load a (double) matrix from a file (usually the coordinates)
double *mesh_load_double (const char *fname, const index cols, index *rows)
{
    FILE *file;
    index cnt, j;
    double a, *data;

    // try open file
    file = fopen(fname,"r");
    // exception handling
    if (file == NULL) return (NULL) ;

    // count entries
    cnt = 0;
    while (fscanf(file,"%lg",&a) != EOF) cnt++;

    // check if matrix is square
    rows[0] = cnt/cols;
    if (  cnt - rows[0] * cols ){
      fclose(file);
      printf("\nmesh_load_double() Error!!! cnt = %g, rows %g\n\n", 
              (double) cnt, (double) rows[0]);
      return (NULL);
    }

    // allocate memory
    data = malloc(cnt * sizeof(double));
    if (!data) return (NULL) ;

    // transmit data to data
    fseek(file,0L,SEEK_SET);
    for (j=0; j<cnt; j++) {
      fscanf(file,"%lg",&(data[j]));
    }

    // close file
    fclose(file);
    return (data);
}


// Load problem from files
// Problem description consists of
// problem.co, problem.el, problem.bd
mesh *mesh_load (char *fname)
{
    // local quantities
    FILE *file;
    char *tmp;
    index cnt, j, a, *data;

    // declare mesh pointer
    mesh *M;
    char buffer[512];

    // initialize zero-sized mesh
    M = mesh_alloc(0,0,0) ;               
    if (!M) return (NULL) ;
    M->nedges = 0;

    // load the coordinates from file
    sprintf(buffer,"%s.co",fname);
    printf("Load coordinates from %s\n",buffer);
    M->coord = mesh_load_double(buffer, 2, &(M->ncoord));

    // load the elements from file
    sprintf(buffer,"%s.el",fname);
    printf("Load elements from %s\n",buffer);
    M->elem = mesh_load_index(buffer, 7, &(M->nelem));

    // load the boundary from file
    sprintf(buffer,"%s.bd",fname);
    printf("Load boundary data from %s\n",buffer);
    M->bdry = mesh_load_index(buffer, 4, &(M->nbdry));

    // exception handling
    return ((!M->coord || !M->elem || !M->bdry) ? mesh_free (M) : M) ;
}
