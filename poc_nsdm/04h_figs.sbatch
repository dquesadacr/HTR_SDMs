#!/bin/bash
#SBATCH --job-name=stage2fig
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --time=00:20:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm-%j.out
POC=/p/projects/isimip/isimip/dannellq/SDMs/poc_nsdm
apptainer exec -B /p -B /tmp -B /home /p/projects/isimip/isimip/dannellq/SDMs/sdm_2026.sif Rscript $POC/07_stage2_figs.R $POC
