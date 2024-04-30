#!/bin/sh

module --force purge
module load singularity/3.8.1

base=`pwd`
export cont="$base"/spt_sdm.sif

cd "$2"/C"$1"
dir=`pwd`

mkdir -p $dir/1_Inputs/2_Predictors/{1_Current/Hist_RCM,2_Projection/RCP85}

cp $base/preds/Past/* $dir/1_Inputs/2_Predictors/1_Current/
cp $base/preds/Hist/* $dir/1_Inputs/2_Predictors/1_Current/Hist_RCM/
cp $base/preds/RCP85/* $dir/1_Inputs/2_Predictors/2_Projection/RCP85/

echo singularity exec -C -B $dir:/Data $cont bash -c "cd /Data;
Rscript colinvar_all.R -f "$1" ; wait
Rscript select_vars_SD.R &
Rscript select_vars_WC_CMIP6.R; wait"

singularity exec -C -B $dir:/Data $cont bash -c "cd /Data;
Rscript colinvar_all.R -f "$1" ; wait
Rscript select_vars_SD.R &
Rscript select_vars_WC_CMIP6.R; wait"

wait

# for j in 2 3 5 10; do # kfolds, set to 10 after testing
for j in 10; do
    cd $base/"$2"
    dir2=$base/"$2"/C"$1"_F"$j"
    cp -a $dir $dir2
    cd $dir2
    for k in alpina montana caryo poa pseudo scabi dianthus meum ranun; do
        echo sbatch -c 16 --qos=medium --mem=58G -J "$k"_C"$1"_F"$j"_"$2" -o $base/logs/"$k"_C"$1"_F"$j"_"$2"_%j.out -e $base/logs/"$k"_C"$1"_F"$j"_"$2"_%j.err --mail-user dannell.quesada@pik-potsdam.de --mail-type END $dir2/run_chain.sh $dir2 $cont $k $j "$3"
        sbatch -c 16 --qos=medium --mem=58G -J "$k"_C"$1"_F"$j"_"$2" -o $base/logs/"$k"_C"$1"_F"$j"_"$2"_%j.out -e $base/logs/"$k"_C"$1"_F"$j"_"$2"_%j.err --mail-user dannell.quesada@pik-potsdam.de --mail-type END $dir2/run_chain.sh $dir2 $cont $k $j "$3"
    done
done
