# HTR_SDMs

This repository contains the source code necessary to reproduce the results of the manuscript "Integrating High-Temporal-Resolution Climate Projections into Species Distribution Models", to be submitted.

## Other resources needed

The following external resources are publicly available to fully reproduce the results:

- The `spt_sdm.sif` apptainer container ([DOI](https://doi.org/10.5281/zenodo.8352758))
  - This container holds all the necessary software and packages to run the code within this repository
- The `preds.zip` preprocessed predictors file ([DOI](https://doi.org/10.5281/zenodo.8404563))

> **Note:** The presented code was run from an HPC environment with SLURM, yet it is possible to adapt to make it run locally.

## Setting everything up and reproducing the results

From the root folder of the project:

1. We download the container and the predictors (and unzip them):

   ```shell
   wget https://zenodo.org/records/8352758/files/spt_sdm.sif?download=1 -O spt_sdm.sif

   wget https://zenodo.org/records/8404563/files/predictors.zip?download=1 -O preds.zip
   unzip preds.zip
   rm preds.zip
   ```

2. We clone this repository, keeping only the necessary files within this folder, and create a `logs` folder:

   ```shell
   git clone --depth 1 https://github.com/dquesadacr/HTR_SDMs.git

   cp HTR_SDMs/train_eval_predict.sh ./
   cp HTR_SDMs/run_colinvar.sh ./
   cp HTR_SDMs/0_Code/brms_SchoD_plots.sh ./

   rm -rf HTR_SDMs
   mkdir -p logs
   ```

3. We submit the jobs to train the SDMs for both the Full version and the spatiotemporal transferability test:

   ```shell
   bash train_eval_predict.sh Full
   wait
   bash train_eval_predict.sh LatLon "_latlon"
   wait
   ```

4. When all the `run_chain.sh` jobs inside the loops of both `train_eval_predict.sh` have been submitted, the postprocessing and plotting jobs are submitted:

   ```shell
   # The following can include unrelated jobs
   afterok=$(squeue --format="%i" -u $USER | tail -n +2 | xargs | tr " " ":")

   sbatch --nodes=1 -c 8 --mem=29GB --qos=medium \
       --dependency=afterok:"$afterok" -J brms_SchoD \
       -o ./logs/brms_SchoD_%j.out -e ./logs/brms_SchoD_%j.err \
       --mail-user dannell.quesada@pik-potsdam.de --mail-type END \
       brms_SchoD_plots.sh "Full/C1_F10" "$cont"
   ```

> [!NOTE]
> The `submit_pipeline.sh` script contains these commands.

## Citation

The bibliography to cite when using either this repository, the _apptainer_ container or the preprocessed predictors will be added in due time.
