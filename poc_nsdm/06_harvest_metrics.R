#!/usr/bin/env Rscript
# 06_harvest_metrics.R — harvest the FULL metric panel (SORENSEN, JACCARD, FPB,
# OR, TSS, AUC, BOYCE/CBI, IMAE + TPR/TNR, mean/sd) that flexsdm already computed
# inside every saved fit object of the paper's C1_F10 run tree, but the paper
# mostly did not report. Read-only over 2_Outputs; writes two CSVs.
# Usage: Rscript 06_harvest_metrics.R <RUN_DIR> <OUT_PREFIX>
suppressMessages(library(dplyr))
args <- commandArgs(TRUE)
RUN <- args[1]; OUTP <- args[2]
files <- list.files(file.path(RUN, "2_Outputs/1_Current"),
                    pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
files <- files[grepl("/(Algorithm|Ensemble)/[^/]+/0_models/", files)]
cat(length(files), "model rds files\n")
perf <- list(); hyper <- list(); nread <- 0
for (f in files) {
  parts <- strsplit(f, "/")[[1]]
  i <- which(parts == "0_models")
  sp <- parts[i + 1]; fit <- sub("[.]rds$", "", parts[i + 2])
  kind <- if (parts[i - 2] == "Algorithm") paste0("alg:", parts[i - 1]) else paste0("ens:", parts[i - 1])
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  nread <- nread + 1
  if (is.null(obj) || is.null(obj$performance)) next
  p <- as.data.frame(obj$performance)
  p$species <- sp; p$kind <- kind; p$run <- fit
  perf[[length(perf) + 1]] <- p
  if (!is.null(obj$hyper_performance)) {
    h <- as.data.frame(obj$hyper_performance)
    h$species <- sp; h$kind <- kind; h$run <- fit
    hyper[[length(hyper) + 1]] <- h
  }
  if (nread %% 50 == 0) { cat(nread, "/", length(files), "\n"); gc(verbose = FALSE) }
}
out <- bind_rows(perf)
write.csv(out, paste0(OUTP, "_performance.csv"), row.names = FALSE)
cat("performance rows:", nrow(out), "->", paste0(OUTP, "_performance.csv\n"))
if (length(hyper)) {
  h <- bind_rows(hyper)
  write.csv(h, paste0(OUTP, "_hyperperformance.csv"), row.names = FALSE)
  cat("hyper rows:", nrow(h), "->", paste0(OUTP, "_hyperperformance.csv\n"))
}
