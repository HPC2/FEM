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
" > $9
}

echo "make mpi"
make mpi > /dev/null
echo "make seq"
make seq > /dev/null

RESULT_NAME="main_square_eq_dof"
HOSTFILE_NAME="hf_main_square_eq_dof"

# 1 1 11; 2 2 10; 4 4 9; 8 8 8
run_experiment 0 1 "./seq_fem 1 1 11 cg $RESULT_NAME"
hostfile 4 0 0 0 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 4 --hostfile $HOSTFILE_NAME mpi_fem 2 2 10 cg $RESULT_NAME"
hostfile 4 4 4 4 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 16 --hostfile $HOSTFILE_NAME mpi_fem 4 4 9 cg $RESULT_NAME"
hostfile 4 4 4 4 4 4 4 4 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 64 --hostfile $HOSTFILE_NAME mpi_fem 8 8 8 cg $RESULT_NAME"

# 1 1 12; 2 2 11; 4 4 10; 8 8 9
run_experiment 0 1 "./seq_fem 1 1 12 cg $RESULT_NAME"
hostfile 4 0 0 0 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 4 --hostfile $HOSTFILE_NAME mpi_fem 2 2 11 cg $RESULT_NAME"
hostfile 4 4 4 4 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 16 --hostfile $HOSTFILE_NAME mpi_fem 4 4 10 cg $RESULT_NAME"
hostfile 4 4 4 4 4 4 4 4 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 64 --hostfile $HOSTFILE_NAME mpi_fem 8 8 9 cg $RESULT_NAME"

# 1 1 13; 2 2 12; 4 4 11; 8 8 10
run_experiment 0 1 "./seq_fem 1 1 13 cg $RESULT_NAME"
hostfile 4 0 0 0 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 4 --hostfile $HOSTFILE_NAME mpi_fem 2 2 12 cg $RESULT_NAME"
hostfile 4 4 4 4 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 16 --hostfile $HOSTFILE_NAME mpi_fem 4 4 11 cg $RESULT_NAME"
hostfile 4 4 4 4 4 4 4 4 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 64 --hostfile $HOSTFILE_NAME mpi_fem 8 8 10 cg $RESULT_NAME"

# 1 1 10; 2 2 9; 4 4 8; 8 8 7
run_experiment 0 1 "./seq_fem 1 1 12 cg $RESULT_NAME"
hostfile 4 0 0 0 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 4 --hostfile $HOSTFILE_NAME mpi_fem 2 2 11 cg $RESULT_NAME"
hostfile 4 4 4 4 0 0 0 0 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 16 --hostfile $HOSTFILE_NAME mpi_fem 4 4 10 cg $RESULT_NAME"
hostfile 4 4 4 4 4 4 4 4 $HOSTFILE_NAME
run_experiment 0 1 "mpirun -np 64 --hostfile $HOSTFILE_NAME mpi_fem 8 8 9 cg $RESULT_NAME"

echo "finished"
exit 0

