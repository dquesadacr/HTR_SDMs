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
preds <- c(cmdArg("p"), cmdArg("preds"), cmdArg("predictors")) # "biovars" # read as BV (bv, biovars) also, and indexes or ind

plot_bool <- c(cmdArg("plot"), cmdArg("pl"))
Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID")) #"centaurea" #
# trimObs2H <- as.logical(c(cmdArg("trimobs"), cmdArg("trim"), cmdArg("TO"))) # trim observations only in during historical period

if(length(plot_bool)==0) {plot_bool <- FALSE}
# if(length(trimObs2H)==0) {trimObs2H <- FALSE}

if(str_detect(preds, "biovars|BV|bv")) {preds <- "biovars"}
if(str_detect(preds, "ind|indexes")) {preds <- "ind"}
if(str_detect(preds, "mix|both")) {preds <- "both"}

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

# READING DATA =================================================================

sapply(c("Hist_RCM", "RCP85"), function(expmt){
  ifelse(expmt == "Hist_RCM", 
         expmt_fold <- "1_Current/Hist_RCM", 
         expmt_fold <- paste0("2_Projection/", expmt))
  
  if(!file.exists(paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/1_con/", 
                         Spec4folders, "/preds_ens_", preds, ".rds"))){
    ens_list <- sapply(list.files(paste0("./2_Outputs/1_Current/Ensemble/meanw/0_models/", Spec4folders),
                                  paste0("^", preds,"*"), full.names = TRUE), simplify = FALSE, readRDS)
    names(ens_list) <- names(ens_list) %>% str_split_i("/", -1) %>% str_remove(".rds")

    preds_list <- sapply(list.files(paste0("./1_Inputs/2_Predictors/", expmt_fold),
                                    paste0("^",preds,"*"), full.names = TRUE) %>%
                           str_subset("colinvar|full|all", negate = TRUE),
                         simplify = FALSE, function(x){
                           pred <- readRDS(x)
                           if (is(pred, "SpatRaster")) {
                             pred <- list(H2000 = terra::unwrap(pred))
                           } else {
                             pred <- lapply(pred, terra::unwrap)
                           }
                           return(pred)
                         })
    names(preds_list) <- names(preds_list) %>% str_split_i("/", -1) %>% str_remove_all(paste0(".rds"))

    # PREDICT ====
    pred_all <- ens_predict(sp_ID = Spec4folders, ens_list, preds_list,
                            plot = plot_bool, expmt_fold = expmt_fold,
                            preds = preds, fileid = preds) %>% unwrap_lists()
  }
})
