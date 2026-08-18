#!/bin/sh

module --force purge

dir=$1
cont=$2
id=$3
folds=$4

echo apptainer exec -B /p -B /home --pwd $dir $cont bash -c "bash full_chain'$5'.sh $id $folds FALSE"

apptainer exec -B /p -B /home --pwd $dir $cont bash -c "bash full_chain'$5'.sh $id $folds FALSE" #TRUE for plotting
