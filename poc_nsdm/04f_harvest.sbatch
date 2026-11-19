#!/bin/bash
#SBATCH --job-name=harvest-m
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%j.out
RUN=/p/projects/isimip/isimip/dannellq/SDMs/Full/C1_F10
POC=/p/projects/isimip/isimip/dannellq/SDMs/poc_nsdm
apptainer exec -B /p -B /tmp -B /home /p/projects/isimip/isimip/dannellq/SDMs/sdm_2026.sif Rscript $POC/06_harvest_metrics.R $RUN $POC/metrics_paper
