#!/bin/bash

sigterm_handler() { 
  echo "Shutdown signal received."
  exit 1
}

## Setup signal trap
trap 'trap " " SIGINT SIGTERM SIGHUP; kill 0; wait; sigterm_handler' SIGINT SIGTERM SIGHUP

function run_experiment() {
    WARMUPS=$1
    RUNS=$2
    EXPERIMENT=$3
    if [ $WARMUPS -gt 0 ]; then
        echo "Warmup experiment '$EXPERIMENT' $WARMUPS times"
        for i in `seq 1 $WARMUPS`; do
            eval "${EXPERIMENT}_warmup"
        done
    fi
    echo "Run experiment '$EXPERIMENT' $RUNS times"
    for i in `seq 1 $RUNS`; do
        eval "${EXPERIMENT}"
    done
}

function hostfile(){
echo "daucher slots=$1
ensinger slots=$2
erhart slots=$3
faulhaber slots=$4
geyer slots=$5
harder slots=$6
krafft slots=$7
magirus slots=$8
" > hostfile
}

echo "make mpi"
make mpi > /dev/null
echo "make seq"
make seq > /dev/null

RESULT_NAME="main_cg_strip"

hostfile 0 0 0 0 0 0 0 2
run_experiment 30 50 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 1 cg $RESULT_NAME"
run_experiment 30 50 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 2 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 3 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 4 cg $RESULT_NAME"
run_experiment 5 10 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 5 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 6 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 7 cg $RESULT_NAME"
run_experiment 0 3 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 8 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 9 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 2 --hostfile hostfile mpi_fem 1 2 10 cg $RESULT_NAME"

hostfile 0 0 0 0 0 0 0 4
run_experiment 30 50 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 1 cg $RESULT_NAME"
run_experiment 30 50 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 2 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 3 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 4 cg $RESULT_NAME"
run_experiment 5 10 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 5 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 6 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 7 cg $RESULT_NAME"
run_experiment 0 3 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 8 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 9 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 4 --hostfile hostfile mpi_fem 1 4 10 cg $RESULT_NAME"

hostfile 0 0 0 0 0 0 4 4
run_experiment 30 50 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 1 cg $RESULT_NAME"
run_experiment 30 50 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 2 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 3 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 4 cg $RESULT_NAME"
run_experiment 5 10 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 5 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 6 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 7 cg $RESULT_NAME"
run_experiment 0 3 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 8 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 9 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 8 --hostfile hostfile mpi_fem 1 8 10 cg $RESULT_NAME"

hostfile 0 0 0 0 4 4 4 4
run_experiment 30 50 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 1 cg $RESULT_NAME"
run_experiment 30 50 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 2 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 3 cg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 4 cg $RESULT_NAME"
run_experiment 5 10 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 5 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 6 cg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 7 cg $RESULT_NAME"
run_experiment 0 3 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 8 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 9 cg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 16 --hostfile hostfile mpi_fem 1 16 10 cg $RESULT_NAME"

exit 0

