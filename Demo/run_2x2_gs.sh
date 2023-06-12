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

RESULT_NAME="main_2x2_gs"

run_experiment 30 50 "mpirun -np 4 mpi_fem 2 2 1 gs $RESULT_NAME"
run_experiment 30 50 "mpirun -np 4 mpi_fem 2 2 2 gs $RESULT_NAME"
run_experiment 10 30 "mpirun -np 4 mpi_fem 2 2 3 gs $RESULT_NAME"
run_experiment 10 30 "mpirun -np 4 mpi_fem 2 2 4 gs $RESULT_NAME"
run_experiment 30 50 "./seq_fem 2 2 1 gs $RESULT_NAME"
run_experiment 30 50 "./seq_fem 2 2 2 gs $RESULT_NAME"
run_experiment 10 30 "./seq_fem 2 2 3 gs $RESULT_NAME"
run_experiment 10 30 "./seq_fem 2 2 4 gs $RESULT_NAME"

exit 0

