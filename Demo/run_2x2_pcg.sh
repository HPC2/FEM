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

echo "make mpi"
make mpi > /dev/null
echo "make seq"
make seq > /dev/null

RESULT_NAME="main_2x2_pcg"

run_experiment 30 50 "mpirun -np 4 mpi_fem 2 2 1 pcg $RESULT_NAME"
run_experiment 30 50 "mpirun -np 4 mpi_fem 2 2 2 pcg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 4 mpi_fem 2 2 3 pcg $RESULT_NAME"
run_experiment 10 30 "mpirun -np 4 mpi_fem 2 2 4 pcg $RESULT_NAME"
run_experiment 5 10 "mpirun -np 4 mpi_fem 2 2 5 pcg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 4 mpi_fem 2 2 6 pcg $RESULT_NAME"
run_experiment 0 5 "mpirun -np 4 mpi_fem 2 2 7 pcg $RESULT_NAME"
run_experiment 0 3 "mpirun -np 4 mpi_fem 2 2 8 pcg $RESULT_NAME"
run_experiment 0 1 "mpirun -np 4 mpi_fem 2 2 9 pcg $RESULT_NAME"
run_experiment 30 50 "./seq_fem 2 2 1 pcg $RESULT_NAME"
run_experiment 30 50 "./seq_fem 2 2 2 pcg $RESULT_NAME"
run_experiment 10 30 "./seq_fem 2 2 3 pcg $RESULT_NAME"
run_experiment 10 30 "./seq_fem 2 2 4 pcg $RESULT_NAME"
run_experiment 5 10 "./seq_fem 2 2 5 pcg $RESULT_NAME"
run_experiment 0 5 "./seq_fem 2 2 6 pcg $RESULT_NAME"
run_experiment 0 5 "./seq_fem 2 2 7 pcg $RESULT_NAME"
run_experiment 0 3 "./seq_fem 2 2 8 pcg $RESULT_NAME"
run_experiment 0 1 "./seq_fem 2 2 9 pcg $RESULT_NAME"

echo "finished"
exit 0

