#!/bin/sh

module --force purge

dir=$1
cont=$2
id=$3
folds=$4
full_chain_suffix=$5

echo apptainer exec -B /p -B /home --pwd $dir $cont bash -c "bash full_chain$full_chain_suffix.sh $id $folds FALSE"

apptainer exec -B /p -B /home --pwd $dir $cont bash -c "bash full_chain$full_chain_suffix.sh $id $folds FALSE" #TRUE for plotting
