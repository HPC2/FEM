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
        echo "Warmup experiment '$EXPERIMENT' $WARMUPS times" >> $LOGFILE
        for i in `seq 1 $WARMUPS`; do
            eval "${EXPERIMENT}_warmup"
        done
    fi
    echo "Run experiment '$EXPERIMENT' $RUNS times" >> $LOGFILE
    for i in `seq 1 $RUNS`; do
        eval "${EXPERIMENT}"
    done
}
function hostfile(){
echo "acker slots="${1}"
daucher slots="${2}"
ensinger slots="${3}"
erhart slots="${4}"
faulhaber slots="${5}"
geyer slots="${6}"
harder slots="${7}"
heim slots="${8}"
krafft slots="${9}"
magirus slots="${10}"
mauch slots="${11}"
multscher slots="${12}"
syrlin slots="${13}"
unseld slots="${14}"
wolbach slots="${15}"
zeitblom slots="${16}"
" > "${17}"
}

DOM="$1"
REF="$2"
DT="$(date +%s)"
RESULT="square_eq_dof_${DOM}_${DOM}_${REF}_$DT"
NP="$(($DOM*$DOM))"
echo $RESULT
LOGFILE="../Result/$RESULT.log"
HOSTFILE="hf_$RESULT"

echo "Logging into $LOGFILE"
echo "Saving into result_$RESULT.csv"

case $NP in
    4) hostfile 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 $HOSTFILE;;
    16) hostfile 0 0 4 4 4 4 0 0 0 0 0 0 0 0 0 0 $HOSTFILE;;
    64) hostfile 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 $HOSTFILE;;
    *) echo "Wrong domain";exit 1;;
esac

run_experiment $3 $4 "mpirun -np 4 --hostfile $HOSTFILE mpi_fem $DOM $DOM $REF cg $RESULT"

echo "finished" >> $LOGFILE
exit 0