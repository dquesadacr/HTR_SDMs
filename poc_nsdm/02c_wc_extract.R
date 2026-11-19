#!/usr/bin/env Rscript
# 02c_wc_extract.R — WorldClim 2.1 base (1970-2000) bioclim, any resolution, cropped to Europe+ box.
# Run inside sdm_2026.sif on HPC (compute node; no internet needed).
# Usage: Rscript 02c_wc_extract.R <unzipped_wc_dir> <out.rds>
suppressMessages(library(terra))
args <- commandArgs(TRUE); DIR <- args[1]; OUT <- args[2]
fs <- list.files(DIR, pattern = "^wc2\\.1_\\S+_bio_\\d+\\.tif$", full.names = TRUE)
stopifnot(length(fs) == 19)
i <- as.integer(sub(".*bio_(\\d+)\\.tif", "\\1", basename(fs)))
r <- rast(fs[i])                      # ordered Bio1..Bio19 (lazy, PackedSpatRaster)
names(r) <- paste0("Bio", 1:19)
r <- crop(r, ext(-15, 28, 32, 64))
saveRDS(r, OUT)
cat("wrote", OUT, "| ext:", as.vector(terra::ext(r)), "| nlyr:", nlyr(r), "| res:", as.vector(res(r)), "\n")
