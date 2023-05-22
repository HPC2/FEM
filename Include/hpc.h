#ifndef _HPC_H
#define _HPC_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include <errno.h>
#include <string.h>

#include <stdbool.h>

#define index ptrdiff_t


// =================== primary HPC routines and data structures ================

typedef struct cs_sparse // matrix in compressed-row/col or triplet form
{
    index nzmax ;     /* maximum number of entries */
    index m ;         /* number of rows */
    index n ;         /* number of columns */
    index *p ;        /* col/row pointers (size n+1) or col indices (size nzmax) */
    index *ind ;      /* row/col indices, size nzmax */
    double *x ;       /* numerical values, size nzmax */
    index nz ;        /* # of entries in triplet matrix, 
                       * -1 for compressed-col, -2 for compressed-row */
} cs ;

typedef struct gem_full /* general matrix form, entries stored */
{
    index m ;         /* number of rows */
    index n ;         /* number of columns */
    index incRow;    /* row increment */
    index incCol;    /* col increment */
    double *x ;       /* numerical values */
} gem ;

typedef struct sed_sparse  /* matrix in sparse matrix in compressed col.    */
{                          /* with extracted diagonal storage form          */
    index nzmax ;     /* maximum number of entries                          */
    index   n ;       /* number of rows/columns                             */
    index  *i ;       /* col pointers and row indices                       */
    double *x ;       /* numerical values, size i[n]                        */
} sed ;

typedef struct mesh_data  /* mesh */
{
    index ncoord ;    /* number of coordinates                             */
    index nelem ;     /* number of elements                                */
    index nedges ;    /* number of edges                                   */
    index nbdry ;     /* number of boundary elements                       */
    index nfixed;     /* number of fixed nodes                             */
    double *coord ;   /* coordinates (x1,y1,x2,y2, ... ,x_ncoord,y_ncoord) */
    index *elem ;     /* elements ([n1,n2,n3,m1,m2,m3,t1], ... )           */
    index *bdry ;     /* bdry ([n1,n2,m1,t1], [n3,n4,m2,t2], ...)          */
    index *edge2no;   /* edge to node ([n1, n2], [n3, n4], ...)            */
    index *fixed;     /* bdry ([n1,n2,m1,t1], [n3,n4,m2,t2], ...)          */
} mesh ;

index* get_local_to_global_numbering(mesh* M, const index rows, 
                                              const index cols,
                                              index refinements);
                                              
typedef struct interface_data
{
      index ncoupl ; /* global number of coupling interfaces */
      index *coupl ; /* global coupling interfaces ([a1,e1,l1,r1,c1], ... ) */
      index *dcoupl ; /* number of nodes on interfaces (no crosspoints) */
      index **icoupl ; /* vector of nodenumbers on interfaces (no crosspts) */
      index *interface_counts; /* num interfaces for each processor*/
      index **interface_ids; /* interface ids for each processor */
} interface_data;

typedef struct coupling_data
{
      index nCrossPts; /* number of all CrossPoints */
      index* crossPts; /* crosspoints */
      index* l2g; /* local to global node numbering */
      index ncoupl ; /* local number of coupling interfaces */
      index *coupl ; /* local coupling interfaces ([a1,e1,l1,r1,c1], ... ) */
      index *dcoupl ; /* number of nodes on interfaces (no crosspoints) */
      index **icoupl ; /* vector of nodenumbers on interfaces (no crosspts) */
} coupling_data;

/* utilities */
void coupling_data_print(coupling_data* coupling, int rank);
void interface_data_write(interface_data *interface_data, char* fname);
coupling_data* mpi_split_interfaces(interface_data* interfaces, index* l2g, int n_nodes);
index* mpi_boundaries(index n_rows, index n_cols, index* global_boundaries);

void *hpc_realloc (void *p, index n, size_t size, index *ok);
double hpc_cumsum (index *p, index *c, index n);
 
sed *sed_alloc (index n, index nzmax, index values);
index sed_realloc (sed *A, index nzmax);
sed *sed_free (sed *A);
sed *sed_done (sed *C, void *w, void *x, index ok);
sed *sed_compress (const cs *A);

gem *sed_to_dense(const sed *A, bool sym);
gem *gem_alloc(index m, index n, index incRow, index incCol);
gem *gem_free(gem *A);

index sed_print (const sed *A, index brief);
index sed_dupl (sed *A);
index sed_gs_constr (const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed, bool forward);
                     
/*
Creates a rectangular mesh suited for dividing into n_rows*n_cols
subprocesses.
@param n_rows Number of rows
@param n_cols Number of cols
@param boundaries 4 x 1 array for boundary condition type in anti-clock wise ordering,
starting with the lower boundary of the rectangle
*/
mesh *mesh_create_rect(index n_rows, index n_cols, index *boundaries);

/*
writes a matrix as file "filename.extension"
@param A Pointer to a matrix
@param rows Number of rows
@param cols Number of columns
@param rowmajor true for row-major matrices, false for col-major
@param name Name of the file
@param extension File extension of the file
*/
void print_matrix(index* A, index rows, index cols, bool rowmajor, char* name, char* extension);

/*
writes the local to global numbering matrix
@param l2g_numbering Pointer to the matrix 
@param n_rows the number of row processors
@param n_cols the number of column processors
@param refinements Number of refinements
*/
void write_l2g(index* l2g_numbering, index n_rows, index n_cols, index refinements);

/*
Writes mesh to .co .el .bd files
@param Mesh mesh to write
@param fname filepath (without file extension)
*/
void mesh_write(mesh *Mesh, char* fname);
interface_data* rect_interface_data(mesh* Mesh, index n_rows, index n_cols, index n_refinments);
mesh *mesh_alloc (index ncoord, index nelem, index nbdry);
mesh *mesh_free (mesh *M);
mesh *mesh_load (char *fname);
index *mesh_getFixed(const index nCoord, const index *bdry, 
                     const index nBdry, index *nFixed);
index mesh_print (const mesh *M, index brief);
mesh *mesh_refine(const mesh *In);
mesh* mesh_multi_refine(mesh* in, int refinements);
index mesh_getEdge2no(const index nElem, const index *Elem, 
                      index *nEdges, index** edge2no);


void stima_laplace(double p1[2], double p2[2], double p3[2],
                   index  typ, double dx[3], double ax[3]);

sed *sed_nz_pattern(mesh *M) ; 
index sed_buildS(mesh *M, sed *T);
void mesh_buildRhs(const mesh *M, double *b, double (*f)(double *, index), 
                   double (*g)(double *, index));

double kappa( double x[2], index typ );
double F_vol( double x[2], index typ );

// ================================ BLAS Functions =============================
// === Level 1 ===
void
dcopy(index n,
      const double *x, index incX,
      double *y, index incY);

void
daxpy(index n, double alpha,
      const double *x, index incX,
      double *y, index incY);

double
ddot(index n,
     const double *x, index incX,
     const double *y, index incY);

index
idamax(index n, const double *x, index incX);

void
dswap(index n, double *x, index incX, double *y, index incY);

void
dscal(index  n,
      double alpha,
      double *x, index incX);

// === Level 2 ===
void
sysed_spmv(double alpha,
           const sed *A,
           const double *x, index incX,
           double beta,
           double *y, index incY);


void
daxpyf(size_t m, double alpha,
       const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
       const double *x, ptrdiff_t incX,
       double *y, ptrdiff_t incY);


void
dgemv_axpyf(size_t m, size_t n,
            double alpha,
            const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
            const double *x, ptrdiff_t incX,
            double *y, ptrdiff_t incY);


void
dgemv_dotf(size_t m, size_t n,
           double alpha,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *x, ptrdiff_t incX,
           double *y, ptrdiff_t incY);


void
dgemv(size_t m, size_t n,
      double alpha,
      const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
      const double *x, ptrdiff_t incX,
      double beta,
      double *y, ptrdiff_t incY);


void
dtrsv(size_t n, bool lower, bool unit,
      const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
      double *x, ptrdiff_t incX);

// === Blas Aux. ===
void
printfDGeMatrix(const char * fmt, size_t m, size_t n,
                const double *A,
                ptrdiff_t incRowA, ptrdiff_t incColA);

void
printDGeMatrix(size_t m, size_t n,
               const double *A,
               ptrdiff_t incRowA, ptrdiff_t incColA);

void
printIGeMatrix(size_t m, size_t n,
               const size_t *A,
               ptrdiff_t incRowA, ptrdiff_t incColA);



#define HPC_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define HPC_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define HPC_CSC(A) (A && (A->nz == -1))
#define HPC_CSR(A) (A && (A->nz == -2))
#define HPC_TRIPLET(A) (A && (A->nz >= 0))
#endif

