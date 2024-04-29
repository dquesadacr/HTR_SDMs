
# Bioclimatic &
Rscript SDM_train.R -f "$2" -a dec -p biovars -id "$1" &
Rscript SDM_train.R -f "$2" -a dec -p biovars -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a dec -p biovars -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a dec -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a dec -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a dec -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a dec -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 1 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 5 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 10 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100
wait

Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 1 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 5 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 10 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -trim &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -trim -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -trim -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -trim -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p biovars -id "$1" -trim -ThinAlg spt -spat_thr 2000 -temp_thr 100
wait

# Indexes &
Rscript SDM_train.R -f "$2" -a dec -p ind -id "$1" &
Rscript SDM_train.R -f "$2" -a dec -p ind -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a dec -p ind -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a dec -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a dec -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a dec -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a dec -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 1 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 5 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 10 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100
wait
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 1 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 5 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 10 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -trim &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -trim -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -trim -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -trim -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p ind -id "$1" -trim -ThinAlg spt -spat_thr 2000 -temp_thr 100
wait

# Mix &
Rscript SDM_train.R -f "$2" -a dec -p both -id "$1" &
Rscript SDM_train.R -f "$2" -a dec -p both -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a dec -p both -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a dec -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a dec -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a dec -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a dec -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 1 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 5 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 10 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100
wait
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 0 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 1 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 5 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 10 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -trim &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -trim -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -trim -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -trim -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p both -id "$1" -trim -ThinAlg spt -spat_thr 2000 -temp_thr 100
wait

# Worldclim &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -trim &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -trim -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -trim -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -trim -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc orig -trim -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -ThinAlg spt -spat_thr 2000 -temp_thr 100 &

Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -trim &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -trim -ThinAlg sp -spat_thr 1000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -trim -ThinAlg sp -spat_thr 2000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -trim -ThinAlg spt -spat_thr 1000 -temp_thr 100 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -trim -ThinAlg spt -spat_thr 2000 -temp_thr 100
wait

# Plot thinning results
Rscript thin_plot.R -id "$1" &

# Predict
Rscript SDM_predict-eval.R -p both -id "$1" -pl "$3" &
Rscript SDM_predict-eval.R -p ind -id "$1" -pl "$3" &
Rscript SDM_predict-eval.R -p biovars -id "$1" -pl "$3" &
Rscript SDM_predict-eval.R -p worldclim -id "$1" -pl "$3" &

# wait

# Project
Rscript SDM_project_WC.R -id "$1" -pl "$3" &
Rscript SDM_project_SD.R -p ind -id "$1" -pl "$3" &
Rscript SDM_project_SD.R -p both -id "$1" -pl "$3" &
Rscript SDM_project_SD.R -p biovars -id "$1" -pl "$3"

wait

# Results plots

Rscript eval_plot.R -id "$1" &

wait
echo

