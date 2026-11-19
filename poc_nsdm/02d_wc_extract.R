#!/usr/bin/env Rscript
# 02d_wc_extract.R — WorldClim 2.1 base (1970-2000) bioclim, any resolution, cropped to Europe+ box.
# Layer-by-layer to disk: keeps peak RAM tiny for the 30s global grids.
# Run inside sdm_2026.sif on HPC. Usage: Rscript 02d_wc_extract.R <unzipped_wc_dir> <out.rds> <crop_dir>
suppressMessages(library(terra))
args <- commandArgs(TRUE); DIR <- args[1]; OUT <- args[2]; CROP <- args[3]
fs <- list.files(DIR, pattern = "^wc2\\.1_\\S+_bio_\\d+\\.tif$", full.names = TRUE)
stopifnot(length(fs) == 19)
i <- as.integer(sub(".*bio_(\\d+)\\.tif", "\\1", basename(fs)))
dir.create(CROP, showWarnings = FALSE, recursive = TRUE)
cf <- character(19)
for (k in 1:19) {
  o <- file.path(CROP, sprintf("bio%02d.tif", k))
  if (!file.exists(o)) {
    r <- crop(rast(fs[i[k]]), ext(-15, 28, 32, 64))
    writeRaster(r, o, gdal = c("COMPRESS=LZW", "TILED=YES"), overwrite = TRUE)
  }
  cf[k] <- o
  cat("bio", k, "done\n"); flush(stdout())
}
r <- rast(cf)
names(r) <- paste0("Bio", 1:19)
saveRDS(r, OUT)
cat("wrote", OUT, "| ext:", as.vector(terra::ext(r)), "| nlyr:", nlyr(r), "| res:", as.vector(res(r)), "\n")
