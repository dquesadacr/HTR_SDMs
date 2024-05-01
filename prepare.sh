#!/bin/sh

module --force purge
module load singularity/3.8.1

wget https://zenodo.org/records/8352758/files/spt_sdm.sif?download=1 -O spt_sdm.sif

dir=`pwd`
export cont="$dir"/spt_sdm.sif

wget https://zenodo.org/records/8404563/files/predictors.zip?download=1 -O preds.zip
unzip preds.zip
rm preds.zip

export repo="HTR_SDMs"
git clone --depth 1 https://github.com/dquesadacr/"$repo".git

cp HTR_SDMs/train_eval_predict.sb ./
cp HTR_SDMs/run_colinvar.sh ./
cp HTR_SDMs/0_Code/brms_SchoD_plots.sh ./

rm -rf HTR_SDMs
mkdir -p logs

bash train_eval_predict.sb Full
wait
bash train_eval_predict.sb LatLon "_latlon"
wait

###

sbatch --nodes=1 -c 16 --mem=58GB --qos=medium -J brms_SchoD -o ./logs/brms_SchoD_%j.out -e ./logs/brms_SchoD_%j.err --mail-user dannell.quesada@tu-dresden.de --mail-type END brms_SchoD_plots.sh "Full/C1_F10" "$cont"

afterok=$(squeue --format="%i" -u $USER | tail -n +2 | xargs | tr " " ":")

sbatch --nodes=1 -c 16 --mem=58GB --qos=medium --dependency=afterok:"$afterok" -J brms_SchoD -o ./logs/brms_SchoD_%j.out -e ./logs/brms_SchoD_%j.err --mail-user dannell.quesada@tu-dresden.de --mail-type END brms_SchoD_plots.sh "Full/C1_F10" "$cont"
