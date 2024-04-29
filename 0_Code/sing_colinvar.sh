

Rscript colinvar_all.R -f "$1"
wait

Rscript select_vars_SD.R &
Rscript select_vars_WC_CMIP6.R

wait
echo
