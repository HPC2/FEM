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

RESULT_NAME="main_v3"

declare -a SOLVER=(
    "jacobi"
    "cg"
    "pcg"
)

for s in "${SOLVER[@]}";do
    if [ $s == "jacobi" ]; then
        warmup_low=0
        warmup_high=5
    else
        warmup_low=5
        warmup_high=15
    fi
    # MPI 4x3
    hostfile 0 0 0 0 0 0 0 12
    run_experiment 30 50 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 0 $s $RESULT_NAME"
    run_experiment 30 50 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 1 $s $RESULT_NAME"
    run_experiment 10 30 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 2 $s $RESULT_NAME"
    run_experiment $warmup_high 15 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 3 $s $RESULT_NAME"
    run_experiment $warmup_low 15 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 4 $s $RESULT_NAME"
    run_experiment $warmup_low 10 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 5 $s $RESULT_NAME"
    run_experiment $warmup_low 5 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 6 $s $RESULT_NAME"
    # MPI 12x1
    hostfile 0 0 0 0 0 0 0 12
    run_experiment 30 50 "mpirun -np 12 --hostfile hostfile mpi_fem 12 1 0 $s $RESULT_NAME"
    run_experiment 30 50 "mpirun -np 12 --hostfile hostfile mpi_fem 12 1 1 $s $RESULT_NAME"
    run_experiment 10 30 "mpirun -np 12 --hostfile hostfile mpi_fem 12 1 2 $s $RESULT_NAME"
    run_experiment $warmup_high 15 "mpirun -np 12 --hostfile hostfile mpi_fem 12 1 3 $s $RESULT_NAME"
    run_experiment $warmup_low 15 "mpirun -np 12 --hostfile hostfile mpi_fem 12 1 4 $s $RESULT_NAME"
    run_experiment $warmup_low 10 "mpirun -np 12 --hostfile hostfile mpi_fem 12 1 5 $s $RESULT_NAME"
    run_experiment $warmup_low 5 "mpirun -np 12 --hostfile hostfile mpi_fem 12 1 6 $s $RESULT_NAME"
    # SEQ 4x3
    run_experiment 30 50 "./seq_fem 4 3 0 $s $RESULT_NAME"
    run_experiment 30 50 "./seq_fem 4 3 1 $s $RESULT_NAME"
    run_experiment $warmup_high 30 "./seq_fem 4 3 2 $s $RESULT_NAME"
    run_experiment $warmup_high 15 "./seq_fem 4 3 3 $s $RESULT_NAME"
    run_experiment $warmup_low 15 "./seq_fem 4 3 4 $s $RESULT_NAME"
    run_experiment $warmup_low 10 "./seq_fem 4 3 5 $s $RESULT_NAME"
    run_experiment $warmup_low 5 "./seq_fem 4 3 6 $s $RESULT_NAME"
    # SEQ 1x1
    run_experiment $warmup_high 20 "./seq_fem 1 1 4 $s $RESULT_NAME"
    run_experiment $warmup_high 20 "./seq_fem 1 1 5 $s $RESULT_NAME"
    run_experiment $warmup_low 20 "./seq_fem 1 1 6 $s $RESULT_NAME"
    run_experiment $warmup_low 5 "./seq_fem 1 1 7 $s $RESULT_NAME"
    # MPI 2x2
    hostfile 0 0 0 0 0 0 0 4
    run_experiment 30 50 "mpirun -np 4 --hostfile hostfile mpi_fem 2 2 3 $s $RESULT_NAME"
    run_experiment $warmup_high 20 "mpirun -np 4 --hostfile hostfile mpi_fem 2 2 4 $s $RESULT_NAME"
    run_experiment $warmup_high 20 "mpirun -np 4 --hostfile hostfile mpi_fem 2 2 5 $s $RESULT_NAME"
    run_experiment $warmup_low 20 "mpirun -np 4 --hostfile hostfile mpi_fem 2 2 6 $s $RESULT_NAME"
    # MPI 4x4
    hostfile 0 0 0 0 0 0 8 8
    run_experiment 30 50 "mpirun -np 16 --hostfile hostfile mpi_fem 4 4 2 $s $RESULT_NAME"
    run_experiment 10 30 "mpirun -np 16 --hostfile hostfile mpi_fem 4 4 3 $s $RESULT_NAME"
    run_experiment $warmup_high 20 "mpirun -np 16 --hostfile hostfile mpi_fem 4 4 4 $s $RESULT_NAME"
    run_experiment $warmup_high 20 "mpirun -np 16 --hostfile hostfile mpi_fem 4 4 5 $s $RESULT_NAME"
    # MPI 8x8
    hostfile 8 8 8 8 8 8 8 8
    run_experiment 30 50 "mpirun -np 64 --hostfile hostfile mpi_fem 8 8 1 $s $RESULT_NAME"
    run_experiment 20 50 "mpirun -np 64 --hostfile hostfile mpi_fem 8 8 2 $s $RESULT_NAME"
    run_experiment $warmup_high 30 "mpirun -np 64 --hostfile hostfile mpi_fem 8 8 3 $s $RESULT_NAME"
    run_experiment $warmup_high 10 "mpirun -np 64 --hostfile hostfile mpi_fem 8 8 4 $s $RESULT_NAME"
done


hostfile 0 0 0 0 0 0 0 12
run_experiment 30 50 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 1 gs $RESULT_NAME"
run_experiment 10 30 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 2 gs $RESULT_NAME"
run_experiment 5 15 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 3 gs $RESULT_NAME"
run_experiment 0 15 "mpirun -np 12 --hostfile hostfile mpi_fem 4 3 4 gs $RESULT_NAME"

run_experiment 30 50 "./seq_fem 4 3 0 gs $RESULT_NAME"
run_experiment 30 50 "./seq_fem 4 3 1 gs $RESULT_NAME"
run_experiment 5 30 "./seq_fem 4 3 2 gs $RESULT_NAME"
run_experiment 5 15 "./seq_fem 4 3 3 gs $RESULT_NAME"
run_experiment 0 15 "./seq_fem 4 3 4 gs $RESULT_NAME"

exit 0






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
        #for refs in 6 7;do
        #    run_experiment 0 5 "./seq_fem $dom $refs $solver $RESULT_NAME"
        #done
    done
done


#---------- MPI PROBLEMS -----
declare -a MPI_DOMAINS=(
    "4 3"
)

for dom in "${MPI_DOMAINS[@]}";do
    N_PROCESSORS="$((${dom// /*}))"
    hostfile 0 0 0 0 0 0 0 $N_PROCESSORS
    for solver in "${MPI_SOLVER[@]}";do
        for refs in 0 1 2 3;do
            run_experiment 10 20 "mpirun -q -np $N_PROCESSORS --hostfile hostfile mpi_fem $dom $refs $solver $RESULT_NAME"
        done
        for refs in 4 5;do
            run_experiment 2 10 "mpirun -np $N_PROCESSORS --hostfile hostfile mpi_fem $dom $refs $solver $RESULT_NAME"
        done
        #for refs in 6 7;do
        #    run_experiment 0 5 "mpirun -np $N_PROCESSORS --hostfile hostfile mpi_fem $dom $refs $solver $RESULT_NAME"
        #done
    done
done