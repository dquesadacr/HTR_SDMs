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

## Load modified functions from flexsdm

source("mods.R")
source("aux_proj.R")

# ARGUMENTS TO PARSE =======
preds <- "wc"
plot_bool <- c(cmdArg("plot"), cmdArg("pl"))
Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID"))

if(length(plot_bool)==0) {plot_bool <- FALSE}

if(str_detect(preds, "biovars|BV|bv")) {preds <- "biovars"}
if(str_detect(preds, "wc|WC|worldclim")) {preds <- "worldclim"}
if(str_detect(preds, "ind|indexes")) {preds <- "ind"}
if(str_detect(preds, "mix|both")) {preds <- "both"}

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

# READING DATA =================================================================

#expmt <- "RCP85"
expmt <- "SSP585" # experiment for WC

ens_list <- ensemble_parse(Spec4folders, preds)

preds_list <- predictors_parse(preds, paste0("2_Projection/", expmt))

# PREDICT ====

pred_all <- ens_predict(sp_ID = Spec4folders, ens_list, preds_list,
                        plot = plot_bool,
                        expmt_fold = paste0("2_Projection/", expmt),
                        preds = preds, fileid = preds) %>%
                        unwrap_lists()
