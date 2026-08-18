#!/bin/sh

cd $1
mkdir -p /p/tmp/"$USER"/sdm_tmp
singularity exec -C -B /p:/p -B `pwd`:`pwd` -B /p/tmp/"$USER"/sdm_tmp:/tmp $2 bash -c "cd `pwd`;

Rscript Sorensen_df_ind.R

Rscript brms_all_plot.R

for i in alpina montana caryo poa pseudo scabi dianthus meum ranun ; do
    Rscript SchoD_proj_fix.R -id "'$i'"
done

Rscript final_proj_plots_filt.R
Rscript SchoD_spatial_trans.R
"
