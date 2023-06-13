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
        echo "Warmup experiment '$EXPERIMENT' $WARMUPS times" > $LOGFILE
        for i in `seq 1 $WARMUPS`; do
            eval "${EXPERIMENT}_warmup"
        done
    fi
    echo "Run experiment '$EXPERIMENT' $RUNS times" > $LOGFILE
    for i in `seq 1 $RUNS`; do
        eval "${EXPERIMENT}"
    done
}

DT="$(date +%s)"
RESULT="square_eq_dof_1_1_$1_$DT"
LOGFILE="../Result/$RESULT.log"

echo "Logging into $LOGFILE"
echo "Saving into result_$RESULT.csv"

#run_experiment $2 $3 "./seq_fem 1 1 $1 cg $RESULT_NAME"
echo "$2 $3 ./seq_fem 1 1 $1 cg $RESULT_NAME"
echo "finished" > $LOGFILE
exit 0

