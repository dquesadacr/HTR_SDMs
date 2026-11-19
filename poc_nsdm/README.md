# poc_nsdm — niche-truncation mitigation PoC (N-SDM-lite)

Stress test behind the niche-truncation reframing: can broad-scale occurrence data
(GBIF) mitigate the small-training-window problem shown for the red-list
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

CBI addendum (2026-11-19, job 3341026, `03b_cbi.R` → `results_cbi_30s.csv`):
Continuous Boyce Index via `flexsdm:::boyce(pres, contrast, n_bins=101,
n_width=0.1)` (predictions at test presences vs test pseudo-absences;
sliding-window P:E ratio, Spearman vs window midpoints). The script recomputes
the split AUCs and they match `results_30s.csv` to 1e-15 (split identity proven).

| species | CBI_A | CBI_B | CBI_C |
|---|---|---|---|
| Arnica montana | 0.965 | 0.847 | **0.977** |
| Cicerbita alpina | 0.881 | 0.441 | 0.805 |
| Meum athamanticum | 0.955 | 0.758 | **0.959** |

CBI is threshold-free, so it isolates continuous ranking from the shared-Youden
penalty that drags B/C in the binary metrics. It sharpens the B verdict
(Cicerbita 0.89 AUC but 0.44 CBI on the Saxony landscape) and mildly favours the
ensemble for the well-sampled species (C ≥ A), while C dilutes toward B for
Cicerbita. Note: flexsdm's `boyce` is a windowed P:E Spearman variant of the
continuous Boyce index (Hirzel et al. 2006, after Boyce et al. 2002); for the
manuscript, implement the published P/E-vs-bin-midpoint form and cross-check
against biomod2's CBI.

## Stage 2 — paper machinery under N-SDM nesting (Guisan strategies 3 + 4)

Tier A = the paper's trained C1_F10 `biovars_H2000_ObsTrim_NoThin` meanw ensembles
(loaded, not refit). Tier B = GBIF whole-range (Europe minus Saxony box) on WorldClim
2.1 H2k 30 arcsec, fit with the paper's verbatim six `fit_*/tune_*` grids + meanw,
10-fold part_random CV, env_const PA (envelope on the 2.5' stack). C3w/C3s = strategy 3
(Sorensen-weighted / simple mean of tier suitabilities); C4 = strategy 4 (B suitability
as covariate, A machinery refit). All evaluated on A's own Ext (2001-2015) points at
max_sens_spec; full panel (7 thresholds x all sdm_eval metrics) in results_stage2/.
Harness validity: A rows reproduce the paper's stored Ext evaluations exactly
(verified for all three species).

max_sens_spec (SOR / TSS / AUC / BOYCE / IMAE):

| species (n)        | A_paper | B_broad | C3w   | C3s   | C4    |
|--------------------|---------|---------|-------|-------|-------|
| Arnica (1405+1405) | .955/.910/.985/.841/.839 | **.962/.925**/.985/.923/.655 | .957/.913/**.987**/**.965** | .957/.913/.987/.964 | .950/.898/.984/.823/.800 |
| Cicerbita (33+33)  | **.918/.849**/.955/.865/**.823** | .914/.818/.914/**.915**/.771 | .914/.818/.945/.910 | .906/.818/.945/.910 | .900/.818/**.961**/.881/.794 |
| Meum (2922+2922)   | **.939/.873/.985**/.983/.821 | .735/.519/.854/.863/.649 | .920/.838/.971/**.993** | .921/.841/.972/**.993** | .934/.861/.980/.975/.776 |

Readings: (1) BOYCE (CBI family) improves over A via the C3 combination in ALL three
species (biggest: Arnica .841 -> .965) — the whole-range tier repairs the continuous
suitability surface of the truncated fine-grain model. (2) In the well-sampled Arnica,
B alone and both combos match or beat A on the binary metrics (SOR/TSS) — a genuine
N-SDM gain on the paper's own test set. (3) In weak-B Meum, C3 gracefully degrades
toward A (0.92 vs 0.94 SOR) while still topping BOYCE — robustness argument.
(4) C4 tops Cicerbita AUC (.961). C3w ~= C3s (CV-Sorensen weights near .5/.5).
Files: stage2_<sp>/ on HPC (fitB_* rds, suitability_surfaces.rds, perfB CV csv).
