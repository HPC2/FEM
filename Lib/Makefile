# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
# resp.  
#       make VERBOSE=1
#       make CC=gcc OPENMP=1
#       make CC=/opt/local/bin/clang-mp-3.7 OPT=1 OPENMP=1
#
# To run a demo using the library 
#       cd ../Demo ; make
#

CC = mpicc

LIBRARY = libhpc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC
I = -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
 sed_buildS.o \
 sed_dupl.o \
 sed_print.o \
 sed_gs_constr.o \
 sed_util.o  \
 gem_util.o  \
 hpc_cumsum.o \
 hpc_malloc.o \
 coupling_util.o \
 interface_util.o \
 mesh_buildRhs.o \
 mesh_load.o  \
 mesh_create_rect.o  \
 mesh_getEdge2no.o \
 mesh_print.o  \
 mesh_refine.o  \
 mesh_getFixed.o  \
 mesh_getNumbering.o \
 mesh_mpi_boundaries.o \
 mesh_util.o \
 blas_level1.o \
 blas_level2.o \
 blas_aux.o \
 mpi_split_interfaces.o \
 mpi_sum.o \
 matrix_print.o \
 numbering_write_l2g.o\
 mpi_dot.o \
 mpi_assemble.o \
 comm_buffer_util.o \
 mpi_jacobi.o \
 mpi_cg.o \
 mpi_pcg.o \
 seq_cg.o \
 seq_pcg.o \
 seq_jacobi.o \
 seq_gs.o \
 mesh_flip_edge.o \
 result.o 


$(HPC): ../Include/hpc.h Makefile

%.o: ../Source/%.c ../Include/hpc.h
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(HPC)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install archive in this directory
install: $(AR_TARGET)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

