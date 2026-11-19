#!/bin/bash
#SBATCH --job-name=stage2
#SBATCH --partition=standard
#SBATCH --qos=medium
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
RUN=/p/projects/isimip/isimip/dannellq/SDMs/Full/C1_F10
POC=/p/projects/isimip/isimip/dannellq/SDMs/poc_nsdm
ENS_RUN="$1"   # e.g. biovars_H2000_ObsTrim_NoThin  (passed via: sbatch --export=NONE 04g_stage2.sbatch ENS)
SPECIES="$2"
apptainer exec -B /p -B /tmp -B /home /p/projects/isimip/isimip/dannellq/SDMs/sdm_2026.sif Rscript $POC/05_stage2.R $RUN $POC "$SPECIES" "$ENS_RUN"
