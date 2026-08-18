
#!/usr/bin/env bash
# ============================================================
# full_chain.sh — Train, predict, and evaluate SDM models
#
# Usage: ./full_chain.sh <spec_id> <folds> <plot>
#   spec_id  : species identifier (string)
#   folds    : number of cross-validation folds (integer)
#   plot     : enable plotting (TRUE/FALSE)
# ============================================================

spec_id="$1"
folds="$2"
plot="$3"

# ============================================================
# Training phase — loops over: section → algo → thinning config
# ============================================================

# --- Bioclimatic, Indexes, Mix (biovars / ind / both) ---
for p in biovars ind both; do
  # Algorithms that apply to all three sections
  for a in dec run10 H2000; do
    # Base call (no thinning)
    Rscript SDM_train.R -f "$folds" -a "$a" -p "$p" -id "$spec_id" &

    # Spatial thinning
    for sp_thr in 1000 2000; do
      Rscript SDM_train.R -f "$folds" -a "$a" -p "$p" -id "$spec_id" \
        -ThinAlg sp -spat_thr "$sp_thr" &
    done

    # Spatial + temporal thinning — temp_thr depends on algorithm
    if [ "$a" = "dec" ]; then
      temps=(0 100)
    elif [ "$a" = "run10" ]; then
      temps=(0 1 5 10 100)
    else   # H2000
      temps=(100)
    fi
    for sp_thr in 1000 2000; do
      for t in "${temps[@]}"; do
        Rscript SDM_train.R -f "$folds" -a "$a" -p "$p" -id "$spec_id" \
          -ThinAlg spt -spat_thr "$sp_thr" -temp_thr "$t" &
      done
    done
  done
  wait

  # H2000 with -trim flag (separate block)
  for a in H2000; do
    Rscript SDM_train.R -f "$folds" -a "$a" -p "$p" -id "$spec_id" -trim &
    for sp_thr in 1000 2000; do
      Rscript SDM_train.R -f "$folds" -a "$a" -p "$p" -id "$spec_id" \
        -trim -ThinAlg sp -spat_thr "$sp_thr" &
      Rscript SDM_train.R -f "$folds" -a "$a" -p "$p" -id "$spec_id" \
        -trim -ThinAlg spt -spat_thr "$sp_thr" -temp_thr 100 &
    done
  done
  wait
done

# --- Worldclim (only H2000, with -wc orig/reproj) ---
for wc in orig reproj; do
  Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
    -wc "$wc" &
  for sp_thr in 1000 2000; do
    Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
      -wc "$wc" -ThinAlg sp -spat_thr "$sp_thr" &
    Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
      -wc "$wc" -ThinAlg spt -spat_thr "$sp_thr" -temp_thr 100 &
  done

  # -trim variant
  Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
    -wc "$wc" -trim &
  for sp_thr in 1000 2000; do
    Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
      -wc "$wc" -trim -ThinAlg sp -spat_thr "$sp_thr" &
    Rscript SDM_train.R -f "$folds" -a H2000 -p worldclim -id "$spec_id" \
      -wc "$wc" -trim -ThinAlg spt -spat_thr "$sp_thr" -temp_thr 100 &
  done
  wait
done
wait

# ============================================================
# Plot thinning results
# ============================================================
Rscript thin_plot.R -id "$spec_id" &

# ============================================================
# Predict phase
# ============================================================
for p in both ind biovars worldclim; do
  Rscript SDM_predict-eval.R -p "$p" -id "$spec_id" -pl "$plot" &
done

# ============================================================
# Project phase
# ============================================================
Rscript SDM_project_WC.R -id "$spec_id" -pl "$plot" &
for p in ind both biovars; do
  Rscript SDM_project_SD.R -p "$p" -id "$spec_id" -pl "$plot" &
done

wait

# ============================================================
# Results plots
# ============================================================
Rscript eval_plot.R -id "$spec_id" &
wait
echo

