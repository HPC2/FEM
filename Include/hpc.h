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
#include <sys/time.h>

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

/*
creates a rowmajor matrix with each row containing the local to global numbering for each processor
and return the pointer to this matrix
@param M The global mesh
@param rows The number of rows
@param cols The number of columns
@param refinements The number of refinements
*/
index* get_local_to_global_numbering(mesh* M, const index rows, 
                                              const index cols,
                                              index refinements);

/*
creates a vector of length n_global_nodes with the global to local numbering for one processor
If the node exists, g2l[flobal_node_number] return the local node number
If the node doesn't exist, g2l[global_node_nr] returns -1
@param local_numbering one row (for one processor) of the local-to-global matrix
@param nof_local_nodes The number of local nodes
@param nof_global_nodes The number of global nodes
*/
index* get_global_to_local_numbering(index* local_numbering, 
                                     index  nof_local_nodes,
                                     index  nof_global_nodes);
                                              
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
      index n_global_nodes;
      index n_global_cp; // number of global crosspoints
      index n_local_cp; // number of local crosspoints
      index* crossPts; /* crosspoints */
      index n_v_nodes; // number of interface nodes
      index* v_nodes; // id's of interface nodes
      index n_i_nodes; // number of inner nodes
      index* i_nodes;  // ids of inner nodes
      index* l2g; /* local to global node numbering */
      index ncoupl ; /* local number of coupling interfaces */
      index *coupl ; /* local coupling interfaces ([a1,e1,l1,r1,c1], ... ) */
      index *coupl_sorted; // Sorted coupl by color. Allways length 4. If color does not exits: [0,0,0,0,-1]
      index *dcoupl ; /* number of nodes on interfaces (no crosspoints) */
      index **icoupl ; /* vector of nodenumbers on interfaces (no crosspts) */
      index **icoupl_sorted ; // Same as concept as coupl_sorted
} coupling_data;

// buffers for mpi communication
typedef struct comm_buffers
{
      index n_global_cp; // number of crosspoints
      double* cp_buffer; // buffer for crosspoints
      double* if_buffer1; // send-buffer for interface nodes
      double* if_buffer2; // recv-buffer for interface nodes
} comm_buffers;

double dot_parallel(double* v_i, double* w_i, index n);
double* mpi_assemble_A (sed* A_loc, coupling_data* coupling);
double* mpi_assemble_t2_vec(coupling_data* coupling, double* local_x, index n);
double* mpi_assemble_t1_vec(coupling_data* coupling, double* local_x, index n);


/* utilities */
void sort_coupl(coupling_data* coupling);
void sort_icoupl(coupling_data* coupling);
void coupling_data_print(coupling_data* coupling, int rank);
void interface_data_write(interface_data *interface_data, char* fname);
coupling_data* mpi_split_interfaces(interface_data* interfaces, index* l2g, int n_nodes, index n_global_nodes);
index* mpi_boundaries(index n_rows, index n_cols, index* global_boundaries);

/*
sums up all node-values of the crosspoints by communicating with everyone
@param coupling A pointer to the coupling_data struct
@param x A pointer to the vector containing all local node values for that processor
@param cp_buffer A buffer for sending the node values and receiving the summed up node values. 
It has size (m+1)(n+1) as the crosspoint nodes have the lowest node numbers
*/
void mpi_sum_crosspoints(coupling_data* coupling, double* x, double* cp_buffer);

/*
sums up all node-values of the interfaces by communicating with all 4 neighbors
@param coupling A pointer to the coupling_data struct
@param x A pointer to the vector containing all local node values for that processor
@parma if_buffer_send A buffer for sending the interface node values
@param if_buffer_recv A buffer for receiving the interface node values
*/
void mpi_sum_interfaces(coupling_data* coupling, double* x, double* if_buffer_send, double* if_buffer_recv);

/*
converts a type2 vector into a type1 vector by communicating all crosspoint values and interface values
@param coupling A pointer to the coupling_data struct
@param x A pointer to the vector containing all local node values for that processor
@param cp_buffer A buffer for sending the node values and receiving the summed up node values. 
It has size (m+1)(n+1) as the crosspoint nodes have the lowest node numbers
@parma if_buffer_send A buffer for sending the interface node values
@param if_buffer_recv A buffer for receiving the interface node values
*/
void mpi_convert_type2_to_type1(coupling_data* coupling, double* x, comm_buffers* buffers);
double* accumulate_inv_diag(coupling_data* coupling, sed* A, comm_buffers* buffers);
double mpi_dotprod(index n, double* x, double* y);

index mpi_jacobi(sed* A, coupling_data* coupling, comm_buffers* buffers, mesh* local_mesh, double* x, double* b);
index mpi_cg(sed* A, coupling_data* coupling, comm_buffers* buffers, mesh* local_mesh, double* x, double* b);
index mpi_pcg(sed* A, coupling_data* coupling, comm_buffers* buffers, mesh* local_mesh, double* x, double* b);
index seq_gs(sed* A, mesh* Mesh, double* x, double* b);
index seq_cg(sed* A, mesh* Mesh, double* x, double* b);
index seq_pcg(sed* A, mesh* Mesh, double* x, double* b);
index seq_jacobi(sed* A, mesh* Mesh, double* x, double* b);

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
mesh *mesh_create_rect(index n_rows, index n_cols, index *boundaries, double offset_x, double offset_y, double size_x, double size_y);

/*
writes a matrix as file "filename.extension"
@param A Pointer to a matrix
@param rows Number of rows
@param cols Number of columns
@param rowmajor true for row-major matrices, false for col-major
@param name Name of the file
@param extension File extension of the file
*/
void print_imatrix(index* A, index rows, index cols, bool rowmajor, char* name, char* extension);
void print_dmatrix(double* A, index rows, index cols, bool rowmajor, char* name, char* extension);

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
void mesh_flip_edge(mesh* m);
comm_buffers* alloc_comm_buffers(index n_global_cp, coupling_data* coupling);
void free_comm_buffers(comm_buffers* buffers);
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


void cg_seq(const sed * A, double* f, double *u, double tol, double *u_D, const mesh *M);
void set_dir(const mesh *M, double *u );
void set_dir_u_D(const mesh *M, double *u , const double* u_D );
// ================================ BLAS Functions =============================
// === Level 1 ===
void
dcopy(index n,
      const double *x, index incX,
      double *y, index incY);

void
indexed_dcopy(index n,
      const double *x, index incX,
      double       *y, index incY,
      const index *indices);

void
dmult(index n,
      const double *x, index incX,
      double *y, index incY);

void
indexed_dmult(index n,
      const double *x, index incX,
      double       *y, index incY,
      const index *indices);

void
daxpy(index n, double alpha,
      const double *x, index incX,
      double *y, index incY);

void
indexed_daxpy(index n, double alpha,
      const double *x, index incX,
      double       *y, index incY,
      const index *indices);

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

void
indexed_dscal(index  n,
      double alpha,
      double *x, index incX,
      const index *indices);

// === Level 2 ===
void
sysed_spmv(double alpha,
           const sed *A,
           const double *x, index incX,
           double beta,
           double *y, index incY);

void
indexed_sysed_spmv(
           index m,
           index n, 
           double alpha,
           const sed *A,
           const double *x, index incX,
           double beta,
           double *y, index incY,
           const index *indices_x,
           const index *indices_y);


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
indexed_dgemv(size_t m, size_t n,
      double alpha,
      const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
      const double *x, ptrdiff_t incX,
      double beta,
      double *y, ptrdiff_t incY,
      index *indices_x,
      index *indices_y);


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

//=== RESULT
void result_write(
    char* result_name,
    int n_processors,
    int n_rows,
    int n_cols,
    int n_refs,
    char* solver,
    int n_iter_max,
    int n_iter_min,
    int n_iter_sum,
    int n_coords,
    index* boundaries,
    int dt_build_S,
    int dt_build_rhs,
    int dt_solve
);
struct timeval tv[100];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

#endif

