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
preds <- c(cmdArg("p"), cmdArg("preds"), cmdArg("predictors"))
plot_bool <- c(cmdArg("plot"), cmdArg("pl"))
Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID"))
ow_bool <- c(cmdArg("ow"), cmdArg("o"), cmdArg("overwrite"))

if(length(plot_bool)==0) {plot_bool <- FALSE}
if(length(ow_bool)==0) {ow_bool <- FALSE}

if(str_detect(preds, "biovars|BV|bv")) {preds <- "biovars"}
if(str_detect(preds, "wc|WC|worldclim")) {preds <- "worldclim"}
if(str_detect(preds, "ind|indexes")) {preds <- "ind"}
if(str_detect(preds, "mix|both")) {preds <- "both"}

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

# READING DATA =================================================================

ens_list <- ensemble_parse(Spec4folders, preds)

preds_list <- predictors_parse(preds, "1_Current")

preds_full <- str_remove_all(names(preds_list),
                             (paste0("_", str_split_i(names(preds_list), "_", -1) %>% 
                                       unique))) %>% unique %>% paste0(., "_")
if(length(preds_full)>1) {preds_full <- paste0(preds_full, collapse = "|")}

# Train -- Validation list
TV_list <- sapply(list.files(paste0("./1_Inputs/1_Occurrences/", Spec4folders),
                             ".rds$", full.names = TRUE), simplify = FALSE, readRDS)

names(TV_list) <- names(TV_list) %>% str_split_i("/", -1) %>% str_remove_all(".rds")

# PREDICT ====

pred_all <- ens_predict(sp_ID = Spec4folders, ens_list, preds_list,
                        plot = plot_bool, expmt_fold = "1_Current",
                        preds = preds, fileid = preds) %>% unwrap_lists()

# ENSENSEMBLE EVALUATION ========

eval_all <- eval_df(sp_ID = Spec4folders, pred_all, TV_list,
                    preds_full = preds_full, overwrite = ow_bool)
