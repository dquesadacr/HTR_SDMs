# poc_nsdm — niche-truncation mitigation PoC (N-SDM-lite)

Stress test behind the niche-truncation reframing: can broad-scale occurrence data
(GCIF-style, GBIF) mitigate the small-training-window problem shown for the red-list
species in the HTR_SDMs paper?

Design (3 species: *Arnica montana*, *Cicerbita alpina*, *Meum athamanticum*):

- **A** — local specialist: LfULG Saxon occurrences (1-km deduped, `Jahr` 1970–2000
  trim) + Saxon pseudo-absences. The paper's training window.
- **B** — broad generalist: GBIF occurrences in a Europe box, Saxony excluded by a
  buffered bbox, + European background. Never sees a Saxon record.
- **C** — nested-lite ensemble: mean probability of A and B.

All evaluated on the SAME held-out Saxon points (80/20 split of presences, pseudo-
absences split in matching proportion, one Youden threshold fitted on A's training
data only). Metrics: AUC (Mann–Whitney), TSS, Sørensen. Predictors: WorldClim 2.1
base bioclim (H2k, 1970–2000), Bio1/4/5/6/12/15, two tiers: 2.5 arcmin
(`results.csv`) and 30 arcsec ≈ 1 km (`results_30s.csv`), matching the paper's grain.

Pipeline (run on HPC, `sdm_2026.sif`; run tree `/p/projects/.../SDMs/poc_nsdm/`):

1. `01_pull_gbif.py` — GBIF pull (login node; comma-range syntax, offset cap 10k,
   years 1970–2000, uncertainty ≤ 20 km, 1-km dedupe) → `gbif_*.csv`
2. `02d_wc_extract.R` — WC2.1 zip (2.5m/30s base bio tiers) cropped layer-by-layer
   to Europe+ box → `wc21_bioclim_H2k[_30s].rds`
3. `03_trunc_test.R` — mgcv::GAM + randomForest per training set, metrics, A/B/C
4. `04b_run_all.sbatch` / `04c_run_30s.sbatch` — slurm wrappers

Finding (2026-11-19, both tiers, consistent): with plain WC bioclim the Saxon-only
model does NOT fail even for the rare Cicerbita (n_sax=167; AUC_A ≥ 0.93, TSS ≥ 0.68).
B cannot be thresholded onto a small region at all (TSS 0.15–0.53), and C ≈ A (never
far below). So the paper's species failures were not raw niche truncation on
bioclim predictors — they are tied to the HTR/MI predictor family and the paper's
harder validations (temporal extrapolation), not to training extent per se. Broad
data is at worst harmless (C), at best insurance invisible at this grain. See
`AIOS/history/2026-11-18_sdm-paper_reviewer_dossier.md` for the reviewer context.
