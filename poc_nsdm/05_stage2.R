#!/usr/bin/env Rscript
# 05_stage2.R — N-SDM stage 2: the paper's own machinery under Guisan's nesting.
#
# Tier A (subrange, fine-grain) = the paper's TRAINED models themselves:
#   C1_F10 1_Inputs (occurrences+PA rds), biovars_H2000.rds (HTR BV predictors),
#   trained ensemble rds from 2_Outputs. Nothing refit; zero divergence by
#   construction (P3 reuse).
# Tier B (whole-range, coarse)  = GBIF Europe-minus-Saxony occurrences (stage-1
#   pulls) on the WorldClim 2.1 H2k 30 arcsec tier, fit with the SAME flexsdm
#   calls/grids as SDM_train.R (max/raf/gam/gbm/glm/net + meanw), 10-fold
#   part_random CV, env_const PA (sampling envelope on the 2.5arcmin stack to
#   keep M-GBM tractable at Europe extent; fitting/eval on the 30s stack).
# Combos (both evaluated on A's own external validation set, Ext = 2001-2015
#   presences + their env_const PAs — the paper's own Ext subset):
#   C3w  strategy 3: Sorensen-weighted mean of A-ens and B-ens suitability
#   C3s  strategy 3: simple mean
#   C4   strategy 4: refit of A machinery with B suitability as extra predictor
# Metrics: every sdm_eval column per model/ensemble/subset (SORENSEN, TSS, AUC,
#   KAPPA?, BOYCE, IMAE, JACCARD, FPB, OR, TPR/TNR) -> tidy CSV, nothing dropped.
#
# Usage: Rscript 05_stage2.R <RUN_DIR> <POC_DIR> <SPECIES> <ENS_RUN>
#   RUN_DIR  = /p/projects/.../SDMs/Full/C1_F10
#   POC_DIR  = /p/projects/.../SDMs/poc_nsdm        (stage-1 outputs inside)
#   SPECIES  = "Cicerbita alpina"
#   ENS_RUN  = fit_name of the A ensemble, e.g. biovars_H2000_ObsTrim_...
suppressMessages({library(dplyr); library(terra); library(sf); library(flexsdm); library(tibble)})
a <- commandArgs(TRUE)
RUN <- a[1]; POC <- a[2]; SPECIES <- a[3]; ENS_RUN <- a[4]
CODE <- RUN                               # mods.R lives at the run-tree root
source(file.path(CODE, "mods.R"))         # cmdArg not needed; PA + thin mods
setwd(RUN)                                # SDM_train.R-style relative paths
set.seed(42)
SP4 <- gsub(" ", "_", SPECIES)
OUT <- paste0(POC, "/stage2_", SP4)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FOLDS <- 10
BPRED <- paste0("Bio", c(1, 4, 5, 6, 12, 15))

## ---- Tier A: paper's trained objects -----------------------------------------
get_r <- function(f) { o <- readRDS(f); if (inherits(o, "PackedSpatRaster")) terra::unwrap(o) else o }
predsA <- get_r(paste0(RUN, "/1_Inputs/2_Predictors/1_Current/biovars_H2000.rds"))
occA_f <- list.files(paste0(RUN, "/1_Inputs/1_Occurrences/", SP4), "[.]rds$", full.names = TRUE)
occ_key <- sub("^[a-zA-Z]+_", "", ENS_RUN)   # obs_name drops the predictors prefix
occA_f <- occA_f[grepl(occ_key, occA_f)][1]
stopifnot(nzchar(occA_f))
occA <- readRDS(occA_f)                    # list(Train, Val, Ext)
ensA_f <- paste0(RUN, "/2_Outputs/1_Current/Ensemble/meanw/0_models/", SP4, "/", ENS_RUN, ".rds")
stopifnot(file.exists(ensA_f))
ensA <- readRDS(ensA_f)
sorf <- function(ens) { p <- as.data.frame(ens$performance); as.numeric(p$SORENSEN_mean[1]) }
sorA <- sorf(ensA)

## ---- Tier B: occurrences & predictors ----------------------------------------
gb <- read.csv(paste0(POC, "/gbif_", SP4, ".csv"))[, c("lon", "lat")]
# Saxony exclusion box identical to stage 1 (built from paper presences 1970-2000)
lf <- st_read(paste0(RUN, "/1_Inputs/1_Occurrences/subset_repro.gpkg"), quiet = TRUE) |> st_transform(4326)
lf <- lf[lf$Art_wiss == SPECIES & !is.na(lf$Jahr) & lf$Jahr >= 1970 & lf$Jahr <= 2000, ]
xy <- st_coordinates(lf)
sax <- data.frame(lon = xy[, 1], lat = xy[, 2])
sax <- sax[!duplicated(round(sax, 2)), ]
bx <- range(sax$lon) + c(-0.2, 0.2); by <- range(sax$lat) + c(-0.2, 0.2)
saxbox <- c(bx[1], by[1], bx[2], by[2])
inSax <- gb$lon > saxbox[1] & gb$lon < saxbox[3] & gb$lat > saxbox[2] & gb$lat < saxbox[4]
gbB <- gb[!inSax, ]
cat("B presences:", nrow(gbB), "(excluded in-box:", sum(inSax), ")\n")
PdB <- data.frame(x = gbB$lon, y = gbB$lat, pr_ab = 1L)

predsB  <- get_r(paste0(POC, "/wc21_bioclim_H2k_30s.rds"))[[BPRED]]
predsBa <- get_r(paste0(POC, "/wc21_bioclim_H2k.rds"))[[BPRED]]   # 2.5' PA envelope
calibB <- terra::vect(terra::ext(predsB))
names(predsB) <- names(predsBa) <- paste0("bio", c(1, 4, 5, 6, 12, 15))
nmsB <- names(predsB)

partB <- flexsdm::part_random(data = PdB, pr_ab = "pr_ab",
                              method = c(method = "kfold", folds = FOLDS))
paB <- lapply(unique(partB$.part), function(j) {
  dd <- partB[partB$.part == j, ]
  pa <- sample_pseudoabs_mod(data = dd[, c("x", "y", "pr_ab")], x = "x", y = "y",
                             n = sum(dd$pr_ab == 1),
                             method = c("env_const", env = predsBa),
                             rlayer = predsBa, calibarea = calibB)
  pa$.part <- j; pa
}) %>% bind_rows()
sppB <- bind_rows(partB, paB)
ex_spp <- sdm_extract(data = sppB, x = "x", y = "y", env_layer = predsB, filter_na = TRUE)

## ---- Tier B: the paper's six algorithms + meanw (verbatim grids) ---------------
tmax <- fit_max(data = ex_spp, response = "pr_ab", predictors = nmsB, partition = ".part",
                thr = "max_sens_spec", clamp = FALSE, classes = "default",
                pred_type = "cloglog", regmult = 1); gc(FALSE)
trf  <- tune_raf(data = ex_spp, response = "pr_ab", predictors = nmsB, partition = ".part",
                 grid = expand.grid(mtry = seq(1, 7, 1)), thr = "max_sens_spec",
                 metric = "SORENSEN", n_cores = 1); gc(FALSE)
tgam <- fit_gam(data = ex_spp, response = "pr_ab", predictors = nmsB, partition = ".part",
                thr = "max_sens_spec"); gc(FALSE)
tgbm <- tune_gbm(data = ex_spp, response = "pr_ab", predictors = nmsB, partition = ".part",
                 grid = expand.grid(n.trees = c(20, 50, 100), shrinkage = c(0.1, 0.5, 1),
                                    n.minobsinnode = c(1, 3, 5, 7, 9)),
                 thr = "max_sens_spec", metric = "SORENSEN", n_cores = 1); gc(FALSE)
tglm <- fit_glm(data = ex_spp, response = "pr_ab", predictors = nmsB, partition = ".part",
                thr = "max_sens_spec", poly = 2, inter_order = 0); gc(FALSE)
tnet <- tune_net(data = ex_spp, response = "pr_ab", predictors = nmsB, partition = ".part",
                 grid = expand.grid(size = c(2, 4, 6, 8, 10),
                                    decay = c(0.001, 0.05, 0.1, 1, 3, 4, 5, 10)),
                 thr = "max_sens_spec", metric = "SORENSEN", n_cores = 1); gc(FALSE)
tsvm <- tryCatch(tune_svm(data = ex_spp, response = "pr_ab", predictors = nmsB, partition = ".part",
                 grid = expand.grid(C = c(2, 4, 8, 16, 20), sigma = c(0.01, 0.1, 0.2, 0.3, 0.4)),
                 thr = "max_sens_spec", metric = "SORENSEN", n_cores = 1), error = function(e) NULL)
ms <- list(max = tmax, raf = trf, gam = tgam, gbm = tgbm, glm = tglm, net = tnet, svm = tsvm)
ms <- ms[lengths(ms) != 0]
for (nm in names(ms)) saveRDS(ms[[nm]], paste0(OUT, "/fitB_", nm, ".rds"))
ensB <- fit_ensemble(models = ms, ens_method = "meanw", thr = NULL,
                     thr_model = "max_sens_spec", metric = "SORENSEN")
saveRDS(ensB, paste0(OUT, "/fitB_ensemble_meanw.rds"))
sorB <- sorf(ensB)
cat("CV Sorensen — A:", sorA, "| B:", sorB, "\n")

## ---- Predict both tiers over the Saxon grid -----------------------------------
sa <- sdm_predict(models = ensA, pred = predsA, thr = "max_sens_spec",
                  con_thr = FALSE, predict_area = NULL)[[1]]$meanw$meanw
# A grid is Germany Zone 4 (metres), B stack is WGS84: project A's extent first
predsB_sax <- terra::crop(predsB, terra::extend(terra::projectExtent(sa, crs(predsB)), 0.2))
sb_full <- sdm_predict(models = ensB, pred = predsB_sax, thr = "max_sens_spec",
                       con_thr = FALSE, predict_area = NULL)[[1]]$meanw$meanw
sb <- terra::mask(terra::resample(sb_full, sa, method = "bilinear"), sa)
wA <- sorA / (sorA + sorB)
c3w <- wA * sa + (1 - wA) * sb
c3s <- (sa + sb) / 2

## ---- Strategy 4: B suitability as covariate in A's machinery -------------------
preds4 <- c(predsA, terra::set.names(sb, "biovarB"))
nms4 <- names(preds4)
exA <- sdm_extract(data = occA$Train, x = "x", y = "y", env_layer = preds4, filter_na = TRUE)
s4max <- fit_max(data = exA, response = "pr_ab", predictors = nms4, partition = ".part",
                 thr = "max_sens_spec", clamp = FALSE, classes = "default",
                 pred_type = "cloglog", regmult = 1); gc(FALSE)
s4rf  <- tune_raf(data = exA, response = "pr_ab", predictors = nms4, partition = ".part",
                  grid = expand.grid(mtry = seq(1, 7, 1)), thr = "max_sens_spec",
                  metric = "SORENSEN", n_cores = 1); gc(FALSE)
s4gm  <- fit_gam(data = exA, response = "pr_ab", predictors = nms4, partition = ".part",
                 thr = "max_sens_spec"); gc(FALSE)
s4lm  <- fit_glm(data = exA, response = "pr_ab", predictors = nms4, partition = ".part",
                 thr = "max_sens_spec", poly = 2, inter_order = 0); gc(FALSE)
m4 <- list(max = s4max, raf = s4rf, gam = s4gm, glm = s4lm)
ens4 <- fit_ensemble(models = m4, ens_method = "meanw", thr = NULL,
                     thr_model = "max_sens_spec", metric = "SORENSEN")
saveRDS(ens4, paste0(OUT, "/fitC4_ensemble_meanw.rds"))
s4 <- sdm_predict(models = ens4, pred = preds4, thr = "max_sens_spec",
                  con_thr = FALSE, predict_area = NULL)[[1]]$meanw$meanw

## ---- Evaluate everything on A's Ext (and Val if present) ----------------------
eval_rast <- function(ras, df, tag) {
  ex <- sdm_extract(data = df, x = "x", y = "y", env_layer = ras, filter_na = TRUE)
  ev <- sdm_eval(p = ex %>% filter(pr_ab == 1) %>% pull("meanw") %>% as.numeric,
                 a = ex %>% filter(pr_ab == 0) %>% pull("meanw") %>% as.numeric)
  d <- as.data.frame(ev); d$model <- tag; d$Subset <- NULL; d
}
rows <- list()
for (sub in c("Ext", "Val")) {
  df <- occA[[sub]]
  if (is.null(df) || nrow(df) == 0 || sum(df$pr_ab == 1) < 3) next
  for (tag in c("A_paper", "B_broad", "C3w", "C3s", "C4")) {
    ras <- switch(tag, A_paper = sa, B_broad = sb, C3w = c3w, C3s = c3s, C4 = s4)
    r <- tryCatch(eval_rast(ras, df, tag), error = function(e) { cat(tag, "eval failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(r)) { r$Subset <- sub; r$species <- SPECIES; rows[[length(rows) + 1]] <- r }
  }
}
if (length(rows)) {
  tab <- bind_rows(rows)
  write.csv(tab, paste0(OUT, "/stage2_eval_ext.csv"), row.names = FALSE)
  print(tab %>% select(model, Subset, everything()) |> as.data.frame())
}
# per-model CV performance tables (all metrics, all algorithms) — nothing dropped
allp <- bind_rows(lapply(ms, function(m) { p <- as.data.frame(m$performance); p$species <- SPECIES; p$kind <- "B_flexsdm"; p }),
                  .id = "algo")
write.csv(allp, paste0(OUT, "/stage2_perfB_cv.csv"), row.names = FALSE)
saveRDS(list(sa = terra::wrap(sa), sb = terra::wrap(sb), c3w = terra::wrap(c3w),
             c3s = terra::wrap(c3s), s4 = terra::wrap(s4), sorA = sorA, sorB = sorB),
        paste0(OUT, "/suitability_surfaces.rds"))
cat("STAGE2-DONE", SP4, "\n")
