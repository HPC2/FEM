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

# --- GLOBALS
RESULT_NAME="first_experiment"

declare -a SEQ_SOLVER=(
    "jacobi"
    "gs"
    "cg"
    "pcg"
)
declare -a MPI_SOLVER=(
    "jacobi"
    "cg"
    "pcg"
)

#---------- MPI PROBLEMS -----
declare -a MPI_DOMAINS=(
    "4 3"
)

for dom in "${MPI_DOMAINS[@]}";do
    N_PROCESSORS="$((${dom// /*}))"
    hostfile 0 0 0 0 0 0 0 $N_PROCESSORS
    for solver in "${MPI_SOLVER[@]}";do
        for refs in 0 1 2 3;do
            run_experiment 10 20 "mpirun -np $N_PROCESSORS --hostfile hostfile mpi_fem $dom $refs $solver $RESULT_NAME"
        done
        for refs in 4 5;do
            run_experiment 2 10 "mpirun -np $N_PROCESSORS --hostfile hostfile mpi_fem $dom $refs $solver $RESULT_NAME"
        done
        for refs in 6 7;do
            run_experiment 0 5 "mpirun -np $N_PROCESSORS --hostfile hostfile mpi_fem $dom $refs $solver $RESULT_NAME"
        done
    done
done

#---------- SEQ PROBLEMS -----
declare -a SEQ_DOMAINS=(
    "4 3"
)

for dom in "${SEQ_DOMAINS[@]}";do
    for solver in "${SEQ_SOLVER[@]}";do
        for refs in 0 1 2 3;do
            run_experiment 10 20 "./seq_fem $dom $refs $solver $RESULT_NAME"
        done
        for refs in 4 5;do
            run_experiment 2 10 "./seq_fem $dom $refs $solver $RESULT_NAME"
        done
        for refs in 6 7;do
            run_experiment 0 5 "./seq_fem $dom $refs $solver $RESULT_NAME"
        done
    done
done
