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

bash train_eval_predict.sb Full
bash train_eval_predict.sb LatLon "_latlon"

###

sbatch --nodes=1 -c 16 --mem=58GB --qos=medium -J brms_SchoD -o ./jobs/brms_SchoD_%j.out -e ./jobs/brms_SchoD_%j.err --mail-user dannell.quesada@tu-dresden.de --mail-type END brms_SchoD_plots.sh "Full/C1_F10" "$cont"
