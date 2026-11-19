#!/bin/bash
#SBATCH --job-name=poc30s
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -eu
POC=/p/projects/isimip/isimip/dannellq/SDMs/poc_nsdm
IMG=/p/projects/isimip/isimip/dannellq/SDMs/sdm_2026.sif
WC30=/p/tmp/dannellq/poc_nsdm/wc2.1_30s_bio
CROP=/p/tmp/dannellq/poc_nsdm/wc30s_crop
GPKG=/p/projects/isimip/isimip/dannellq/SDMs/Full/C1_F10/1_Inputs/1_Occurrences/subset_repro.gpkg
# 30 arcsec (~1 km, matches paper) — primary tier
apptainer exec -B /p -B /tmp -B /home $IMG Rscript $POC/02d_wc_extract.R $WC30 $POC/wc21_bioclim_H2k_30s.rds $CROP
apptainer exec -B /p -B /tmp -B /home $IMG Rscript $POC/03_trunc_test.R $POC $POC/wc21_bioclim_H2k_30s.rds $GPKG $POC/results_30s
