#!/bin/bash
#SBATCH --job-name=poc-cbi
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
cd /p/projects/isimip/isimip/dannellq/SDMs/poc_nsdm
apptainer exec -B /p -B /tmp -B /home /p/projects/isimip/isimip/dannellq/SDMs/sdm_2026.sif Rscript 03b_cbi.R . wc21_bioclim_H2k_30s.rds /p/projects/isimip/isimip/dannellq/SDMs/Full/C1_F10/1_Inputs/1_Occurrences/subset_repro.gpkg results_cbi_30s
