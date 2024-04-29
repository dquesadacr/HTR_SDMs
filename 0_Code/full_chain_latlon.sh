
##
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 -trimlon 4615000 -trimlat 5625000 &
Rscript SDM_train.R -f "$2" -a run10 -p biovars -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 5 -trimlon 4615000 -trimlat 5625000 &

Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 -trimlon 4615000 -trimlat 5625000 &
Rscript SDM_train.R -f "$2" -a run10 -p ind -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 5 -trimlon 4615000 -trimlat 5625000 &

Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 0 -trimlon 4615000 -trimlat 5625000 &
Rscript SDM_train.R -f "$2" -a run10 -p both -id "$1" -ThinAlg spt -spat_thr 1000 -temp_thr 5 -trimlon 4615000 -trimlat 5625000 &

Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -ThinAlg sp -spat_thr 1000 -trimlon 4615000 -trimlat 5625000 &
Rscript SDM_train.R -f "$2" -a H2000 -p worldclim -id "$1" -wc reproj -trim -ThinAlg sp -spat_thr 1000 -trimlon 4615000 -trimlat 5625000 &

wait

# Predict
Rscript SDM_predict-eval.R -p both -id "$1" -pl "$3" &
Rscript SDM_predict-eval.R -p ind -id "$1" -pl "$3" &
Rscript SDM_predict-eval.R -p biovars -id "$1" -pl "$3" &
Rscript SDM_predict-eval.R -p worldclim -id "$1" -pl "$3"

wait

# Project
Rscript SDM_project_WC.R -id "$1" -pl "$3" &
Rscript SDM_project_SD.R -p ind -id "$1" -pl "$3" &
Rscript SDM_project_SD.R -p both -id "$1" -pl "$3" &
Rscript SDM_project_SD.R -p biovars -id "$1" -pl "$3"

wait
echo

