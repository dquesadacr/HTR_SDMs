#!/usr/bin/env Rscript
# 03_trunc_test.R — niche-truncation mitigation PoC (N-SDM-lite, 3 species).
# A = Saxony-only (LfULG, truncated niche); B = broad GBIF (Saxony excluded);
# C = nested-lite mean ensemble A+B. All on CHELSA 1979-2008 bioclim at 5 arcmin, 6 predictors.
# Run inside sdm_2026.sif. Usage: Rscript 03_trunc_test.R <poc_dir> <chelesa_box_rds> <gpkg> <out_prefix>
suppressMessages({library(terra); library(sf); library(mgcv); library(randomForest)})
args <- commandArgs(TRUE)
POC <- args[1]; RDS <- args[2]; GPKG <- args[3]; OUTP <- args[4]
set.seed(42)
SP <- c("Arnica montana", "Cicerbita alpina", "Meum athamanticum")
PRED <- paste0("Bio", c(1, 4, 5, 6, 12, 15))
bio <- readRDS(RDS)

lf <- st_read(GPKG, quiet = TRUE) |> st_transform(4326)
lf$sp <- as.character(lf$Art_wiss)
lf <- lf[!is.na(lf$Jahr) & lf$Jahr >= 1970 & lf$Jahr <= 2000, ]   # H2k trimming (paper-consistent)

auc <- function(p, y) {  # Mann-Whitney
  r <- rank(c(p[y == 1], p[y == 0])); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
metrics <- function(p, y, thr) {
  prd <- as.integer(p >= thr)
  tp <- sum(prd == 1 & y == 1); fn <- sum(prd == 0 & y == 1)
  fp <- sum(prd == 1 & y == 0); tn <- sum(prd == 0 & y == 0)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  c(AUC = auc(p, y), TSS = sens + spec - 1,
    SOR = 2 * tp / (fn + 2 * tp + fp))
}
youden <- function(p, y) {
  th <- sort(unique(p)); best <- c(-1, 0.5)
  for (t in quantile(p, seq(0.02, 0.98, 0.02))) {
    prd <- p >= t; tp <- sum(prd & y == 1); fn <- sum(!prd & y == 1)
    fp <- sum(prd & y == 0); tn <- sum(!prd & y == 0)
    j <- tp / (tp + fn) + tn / (tn + fp) - 1
    if (j > best[1]) best <- c(j, t)
  }
  best[2]
}

fit_pred <- function(tr, te) {
  f <- as.formula(paste("pr_ab ~", paste(PRED, collapse = " + ")))
  g <- gam(f, data = tr, family = binomial)
  pf <- as.formula(paste("factor(pr_ab) ~", paste(PRED, collapse = " + ")))
  rf <- randomForest(pf, data = tr, ntree = 500, mtry = min(3, length(PRED)))
  ptr <- (predict(g, tr, type = "response") + predict(rf, tr, type = "prob")[, 2]) / 2
  pte <- (predict(g, te, type = "response") + predict(rf, te, type = "prob")[, 2]) / 2
  list(train = ptr, test = pte)
}

samp <- function(ext_box, n) data.frame(lon = runif(n, ext_box[1], ext_box[3]),
                                        lat = runif(n, ext_box[2], ext_box[4]))
extract_pred <- function(df) {
  v <- vect(df, geom = c("lon", "lat"))
  crs(v) <- crs(bio)
  x <- as.data.frame(terra::extract(bio, v))
  cbind(df, x[, PRED, drop = FALSE])
}

res <- list()
for (sp in SP) {
  gb <- read.csv(sprintf("%s/gbif_%s.csv", POC, gsub(" ", "_", sp)))[, c("lon", "lat")]
  s <- lf[lf$sp == sp, ]
  xy <- st_coordinates(s)
  sax <- data.frame(lon = xy[, 1], lat = xy[, 2])
  sax <- sax[!duplicated(round(sax, 2)), ]                       # 1-km dedupe, like GBIF side
  if (nrow(sax) > 2000) sax <- sax[sample(nrow(sax), 2000), ]
  bx  <- range(sax$lon) + c(-0.2, 0.2); by <- range(sax$lat) + c(-0.2, 0.2)
  saxbox <- c(bx[1], by[1], bx[2], by[2])
  inSax <- gb$lon > saxbox[1] & gb$lon < saxbox[3] & gb$lat > saxbox[2] & gb$lat < saxbox[4]
  gbb <- gb[!inSax, ]
  if (nrow(gbb) > 2000) gbb <- gbb[sample(nrow(gbb), 2000), ]
  # occurrences + predictors
  sP <- extract_pred(sax); sP <- sP[complete.cases(sP[, PRED]), ]
  gP <- extract_pred(gbb); gP <- gP[complete.cases(gP[, PRED]), ]
  nsp <- nrow(sP); if (nsp < 40) { cat(sp, ": too few Saxon points\n"); next }
  itr <- sample(nsp, floor(0.8 * nsp))
  # Saxon pseudoabsences (for A): sample from sax box, split in train/test halves (no PA leakage)
  paA <- extract_pred(samp(saxbox, 6 * nsp)); paA <- paA[complete.cases(paA[, PRED]), ]
  npa <- min(3 * nsp, nrow(paA)); paA <- paA[1:npa, ]
  ipa <- sample(npa, floor(0.8 * npa))
  trA <- rbind(cbind(sP[itr, ], pr_ab = 1), cbind(paA[ipa, ], pr_ab = 0))
  teA <- rbind(cbind(sP[-itr, ], pr_ab = 1), cbind(paA[-ipa, ], pr_ab = 0))
  # broad background (for B)
  bgB <- extract_pred(samp(c(-12, 34, 25, 62), min(4000, 3 * nrow(gP))))
  bgB <- bgB[complete.cases(bgB[, PRED]), ]
  trB <- rbind(cbind(gP, pr_ab = 1), cbind(bgB, pr_ab = 0))
  te <- teA
  m <- function(d) { d <- d[, c(PRED, "pr_ab")]; d$pr_ab <- as.integer(d$pr_ab); d }
  a <- fit_pred(m(trA), m(te)); b <- fit_pred(m(trB), m(te))
  thr <- youden(a$train, trA$pr_ab)
  C <- (a$test + b$test) / 2
  mA <- metrics(a$test, te$pr_ab, thr); mB <- metrics(b$test, te$pr_ab, thr)
  mC <- metrics(C, te$pr_ab, thr)
  names(mA) <- paste0(names(mA), "_A"); names(mB) <- paste0(names(mB), "_B"); names(mC) <- paste0(names(mC), "_C")
  res[[sp]] <- c(species = sp, n_sax = nsp, n_gbif = nrow(gP), n_bgA = nrow(paA), n_bgB = nrow(bgB),
                 thr = thr, mA, mB, mC)
  cat(sp, "A:", paste(round(mA, 3), collapse = "/"), "B:", paste(round(mB, 3), collapse = "/"),
      "C:", paste(round(mC, 3), collapse = "/"), "\n"); flush(stdout())
}
df <- do.call(rbind, res)
write.csv(df, paste0(OUTP, ".csv"))
md <- c("# Niche-truncation PoC results (CHELSA 1979-2008, Bio1/4/5/6/12/15, 80/20 hold-out)",
        sprintf("- generated: %s", Sys.time()),
        "", "| species | n_sax | n_gbif | AUC_A | TSS_A | SOR_A | AUC_B | TSS_B | SOR_B | AUC_C | TSS_C | SOR_C |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
        apply(df, 1, function(r) sprintf("| %s | %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |",
          r["species"], r["n_sax"], r["n_gbif"], as.numeric(r["AUC_A"]), as.numeric(r["TSS_A"]), as.numeric(r["SOR_A"]),
          as.numeric(r["AUC_B"]), as.numeric(r["TSS_B"]), as.numeric(r["SOR_B"]), as.numeric(r["AUC_C"]), as.numeric(r["TSS_C"]), as.numeric(r["SOR_C"]))))
writeLines(md, paste0(OUTP, ".md"))
cat("wrote", paste0(OUTP, ".md\n"))
