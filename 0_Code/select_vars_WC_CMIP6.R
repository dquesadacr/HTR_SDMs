# PACKAGES =====================================================================
## CRAN PACKAGES ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "raster", # moved terra above to avoid masking of several functions
  "tidyverse",
  "terra",
  "sf",
  "tidyterra",
  "ggsci",
  "ggpubr",
  "magrittr",
  "ggnewscale",
  "R.utils"
)
sapply(package_vec, install.load.package)

## NON-CRAN PACKAGES ----
if("flexsdm" %in% rownames(installed.packages()) == FALSE){ # flexsdm check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("sjevelazco/flexsdm")
}
library(flexsdm)

# READING DATA =================================================================

## Climatic Data ----
set.seed(42)

preds2000 <- readRDS("./1_Inputs/2_Predictors/1_Current/worldclim_reproj_H2000.rds")

predsP20 <- sapply(list.files("./1_Inputs/2_Predictors/2_Projection/SSP585", full.names = TRUE, pattern = ".rds"), simplify = FALSE, FUN = function(x){
  readRDS(x) %>% .[[2]] %>% terra::unwrap() %>% select(names(preds2000)) %>% terra::wrap()})

names(predsP20) <- names(predsP20) %>% str_split_i("-", -1) %>% str_remove_all(".rds") %>% as.numeric()

saveRDS(predsP20, paste0("./1_Inputs/2_Predictors/2_Projection/SSP585/worldclim_",
                         list.files("./1_Inputs/2_Predictors/2_Projection/SSP585/",
                                    pattern = "^MPI") %>% 
                           str_split_i(pattern = "_20", 1) %>% unique,
                         ".rds"))
