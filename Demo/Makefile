
CC = gcc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

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

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo0 \
    *.a *.dSYM *.obj *.dll
