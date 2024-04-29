
# echo singularity exec -C -B ./:/Data/ -B /beegfs/ws/1/s9941460-SDM/tmp/:/tmp/ /beegfs/ws/1/s9941460-SDM/sdm.sif3 bash -c "cd /Data; Rscript $2"
cd $1
mkdir -p /p/tmp/"$USER"/sdm_tmp
singularity exec -C -B /p:/p -B `pwd`:`pwd` -B /p/tmp/"$USER"/sdm_tmp:/tmp $2 bash -c "cd `pwd`;

Rscript Sorensen_df_ind.R

Rscript brms_H2k_plot.R

Rscript brms_all_ThTr_plot.R

Rscript brms_TempTrans_plot.R

Rscript brms_R10_plot.R

for i in alpina montana caryo poa pseudo scabi dianthus orchis meum ranun ; do
    Rscript SchoD_proj_fix.R -id $i
done

Rscript final_proj_plots_filt.R
Rscript SchoD_spatial_trans.R
"
