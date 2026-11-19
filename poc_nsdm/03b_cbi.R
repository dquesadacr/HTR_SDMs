#!/usr/bin/env Rscript
# 03b_cbi.R — adds Continuous Boyce Index to the stage-1 A/B/C test.
# Identical to 03_trunc_test.R through all RNG consumption (seed 42), so
# AUC/TSS/SOR must reproduce results_30s.csv exactly — that equality is the
# self-check that the split is unchanged. CBI via flexsdm:::boyce(type="cont"),
# background = TEST pseudo-absences (evaluation PA set).
suppressMessages({library(terra); library(sf); library(mgcv); library(randomForest)})
args <- commandArgs(TRUE)
POC <- args[1]; RDS <- args[2]; GPKG <- args[3]; OUTP <- args[4]
set.seed(42)
SP <- c("Arnica montana", "Cicerbita alpina", "Meum athamanticum")
PRED <- paste0("Bio", c(1, 4, 5, 6, 12, 15))
bio <- readRDS(RDS)
cbi_fun <- getFromNamespace("boyce", "flexsdm")

lf <- st_read(GPKG, quiet = TRUE) |> st_transform(4326)
lf$sp <- as.character(lf$Art_wiss)
lf <- lf[!is.na(lf$Jahr) & lf$Jahr >= 1970 & lf$Jahr <= 2000, ]

cbi <- function(p, te) {
  v <- vect(te, geom = c("lon", "lat")); crs(v) <- crs(bio)
  pres <- v[te$pr_ab == 1, ]; bg <- v[te$pr_ab == 0, ]
  cbi_fun(ecm = p[te$pr_ab == 1], pres = pres, background = bg, type = "cont")
}

fit_pred <- function(tr, te) {
  f <- as.formula(paste("pr_ab ~", paste(PRED, collapse = " + ")))
  g <- gam(f, data = tr, family = binomial)
  pf <- as.formula(paste("factor(pr_ab) ~", paste(PRED, collapse = " + ")))
  rf <- randomForest(pf, data = tr, ntree = 500, mtry = min(3, length(PRED)))
  list(train = (predict(g, tr, type = "response") + predict(rf, tr, type = "prob")[, 2]) / 2,
       test  = (predict(g, te, type = "response") + predict(rf, te, type = "prob")[, 2]) / 2)
}
samp <- function(ext_box, n) data.frame(lon = runif(n, ext_box[1], ext_box[3]),
                                        lat = runif(n, ext_box[2], ext_box[4]))
extract_pred <- function(df) {
  v <- vect(df, geom = c("lon", "lat")); crs(v) <- crs(bio)
  cbind(df, as.data.frame(terra::extract(bio, v))[, PRED, drop = FALSE])
}

res <- list()
for (sp in SP) {
  gb <- read.csv(sprintf("%s/gbif_%s.csv", POC, gsub(" ", "_", sp)))[, c("lon", "lat")]
  s <- lf[lf$sp == sp, ]; xy <- st_coordinates(s)
  sax <- data.frame(lon = xy[, 1], lat = xy[, 2])
  sax <- sax[!duplicated(round(sax, 2)), ]
  if (nrow(sax) > 2000) sax <- sax[sample(nrow(sax), 2000), ]
  bx <- range(sax$lon) + c(-0.2, 0.2); by <- range(sax$lat) + c(-0.2, 0.2)
  saxbox <- c(bx[1], by[1], bx[2], by[2])
  inSax <- gb$lon > saxbox[1] & gb$lon < saxbox[3] & gb$lat > saxbox[2] & gb$lat < saxbox[4]
  gbb <- gb[!inSax, ]
  if (nrow(gbb) > 2000) gbb <- gbb[sample(nrow(gbb), 2000), ]
  sP <- extract_pred(sax); sP <- sP[complete.cases(sP[, PRED]), ]
  gP <- extract_pred(gbb); gP <- gP[complete.cases(gP[, PRED]), ]
  nsp <- nrow(sP); if (nsp < 40) { cat(sp, ": too few Saxon points\n"); next }
  itr <- sample(nsp, floor(0.8 * nsp))
  paA <- extract_pred(samp(saxbox, 6 * nsp)); paA <- paA[complete.cases(paA[, PRED]), ]
  npa <- min(3 * nsp, nrow(paA)); paA <- paA[1:npa, ]
  ipa <- sample(npa, floor(0.8 * npa))
  trA <- rbind(cbind(sP[itr, ], pr_ab = 1), cbind(paA[ipa, ], pr_ab = 0))
  teA <- rbind(cbind(sP[-itr, ], pr_ab = 1), cbind(paA[-ipa, ], pr_ab = 0))
  bgB <- extract_pred(samp(c(-12, 34, 25, 62), min(4000, 3 * nrow(gP))))
  bgB <- bgB[complete.cases(bgB[, PRED]), ]
  trB <- rbind(cbind(gP, pr_ab = 1), cbind(bgB, pr_ab = 0))
  te <- teA
  m <- function(d) { d <- d[, c(PRED, "pr_ab")]; d$pr_ab <- as.integer(d$pr_ab); d }
  a <- fit_pred(m(trA), m(te)); b <- fit_pred(m(trB), m(te))
  C <- (a$test + b$test) / 2
  vals <- c(CBI_A = as.numeric(cbi(a$test, te)),
            CBI_B = as.numeric(cbi(b$test, te)),
            CBI_C = as.numeric(cbi(C, te)))
  res[[sp]] <- c(species = sp, n_sax = nsp, n_gbif = nrow(gP), vals)
  cat(sp, "CBI:", paste(sprintf("%.3f", vals), collapse = " / "), "\n"); flush(stdout())
}
df <- do.call(rbind, res)
write.csv(df, paste0(OUTP, ".csv"))
cat("wrote", paste0(OUTP, ".csv\n"))
