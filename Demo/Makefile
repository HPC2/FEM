
CC = gcc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic


CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib hpc_demo0
	- ./hpc_demo0 problem1 4
	- ./hpc_demo0 problem2 4


lib:
	( cd ../Lib ; $(MAKE) )

hpc_demo0: lib hpc_demo0.c Makefile
	$(CC) $(CF) $(I) -o hpc_demo0 hpc_demo0.c $(CS)

seq: lib seq_fem.c Makefile
	$(CC) $(CF) $(I) -o seq_fem seq_fem.c $(CS)

mpi: lib mpi_fem.c Makefile
	mpicc $(CF) $(I) -o mpi_fem mpi_fem.c $(CS)

cg: lib test_cg_seq.c Makefile
	$(CC) $(CF) $(I) -o cg test_cg_seq.c $(CS)

run_test_cases:
	- mpirun -np 1 --oversubscribe mpi_mesh_split 1 1 1
	- mpirun -np 1 --oversubscribe mpi_mesh_split 1 1 3
	- mpirun -np 4 --oversubscribe mpi_mesh_split 2 2 1
	- mpirun -np 4 --oversubscribe mpi_mesh_split 2 2 3
	- mpirun -np 6 --oversubscribe mpi_mesh_split 2 3 1
	- mpirun -np 6 --oversubscribe mpi_mesh_split 2 3 3
	- mpirun -np 9 --oversubscribe mpi_mesh_split 3 3 1
	- mpirun -np 9 --oversubscribe mpi_mesh_split 3 3 3


clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo0 \
    *.a *.dSYM *.obj *.dll
