#!/bin/bash
#SBATCH --job-name=poc_nsdm
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -eu
POC=/p/projects/isimip/isimip/dannellq/SDMs/poc_nsdm
IMG=/p/projects/isimip/isimip/dannellq/SDMs/sdm_2026.sif
WC=/p/tmp/dannellq/poc_nsdm/wc2.1_2.5_bio
GPKG=/p/projects/isimip/isimip/dannellq/SDMs/Full/C1_F10/1_Inputs/1_Occurrences/subset_repro.gpkg
apptainer exec -B /p -B /tmp -B /home $IMG Rscript $POC/02c_wc_extract.R $WC $POC/wc21_bioclim_H2k.rds
apptainer exec -B /p -B /tmp -B /home $IMG Rscript $POC/03_trunc_test.R $POC $POC/wc21_bioclim_H2k.rds $GPKG $POC/results
