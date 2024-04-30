#!/bin/sh

module --force purge
module load singularity/3.8.1

dir=$1
cont=$2
id=$3
folds=$4

echo singularity exec -C -B $dir:/Data $cont bash -c "cd /Data; bash full_chain'$5'.sh $id $folds FALSE"

singularity exec -C -B $dir:/Data $cont bash -c "cd /Data; bash full_chain'$5'.sh $id $folds FALSE" #TRUE for plotting
