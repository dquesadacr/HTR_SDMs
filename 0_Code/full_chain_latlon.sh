
#!/usr/bin/env bash
# ============================================================
# full_chain_latlon.sh — Lat/lon-trimmed SDM training, predict, evaluate
#
# Usage: ./full_chain_latlon.sh <spec_id> <folds> <plot>
#   spec_id  : species identifier (string)
#   folds    : number of cross-validation folds (integer)
#   plot     : enable plotting (TRUE/FALSE)
# ============================================================

spec_id="$1"
folds="$2"
plot="$3"

TRIMLON=4615000
TRIMLAT=5625000
spat_thr=1000

# ============================================================
# Training phase
# ============================================================

# run10 + spt + spat_thr 1000 + temp_thr 0/5 for biovars, ind, both
for p in biovars ind both; do
  for t in 0 5; do
    Rscript SDM_train.R -f "$folds" -a run10 -p "$p" -id "$spec_id" \
      -ThinAlg spt -spat_thr "$spat_thr" -temp_thr "$t" \
      -trimlon "$TRIMLON" -trimlat "$TRIMLAT" &
  done
done

# H2000 + worldclim + reproj + sp + spat_thr 1000
Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
  -wc reproj -ThinAlg sp -spat_thr "$spat_thr" \
  -trimlon "$TRIMLON" -trimlat "$TRIMLAT" &

# H2000 + worldclim + reproj + trim + sp + spat_thr 1000
Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
  -wc reproj -trim -ThinAlg sp -spat_thr "$spat_thr" \
  -trimlon "$TRIMLON" -trimlat "$TRIMLAT" &

wait

# ============================================================
# Predict phase
# ============================================================
for p in both ind biovars worldclim; do
  Rscript SDM_predict-eval.R -p "$p" -id "$spec_id" -pl "$plot" &
done
wait

# ============================================================
# Project phase
# ============================================================
Rscript SDM_project_WC.R -id "$spec_id" -pl "$plot" &
for p in ind both biovars; do
  Rscript SDM_project_SD.R -p "$p" -id "$spec_id" -pl "$plot" &
done
wait

echo

