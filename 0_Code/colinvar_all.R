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

# ARGUMENTS TO PARSE =======

approach <- "run10"
filt <- as.character(c(cmdArg("filter"), cmdArg("f")))

if(filt=="1"){filt<-" "}
if(filt=="2"){filt<-"Rn_"}
if(filt=="3"){filt<-"days|growing|FFP"}
if(filt=="4"){filt<-"days|growing|FFP|Rn_"}

# DIRECTORIES ==================================================================
Dir.Base <- getwd()
setwd(Dir.Base)

project_directory <- sdm_directory(
  main_dir = Dir.Base,
  projections = c("RCP85"),
  calibration_area = TRUE,
  algorithm = c("all"),
  ensemble = c("meanw"),
  threshold = TRUE,
  return_vector = TRUE
)

# READING DATA =================================================================

## Climatic Data ----
set.seed(42)

# Bioclimatic variables derived from ReKIS ====
preds2000_bv <- readRDS("./1_Inputs/2_Predictors/1_Current/biovars_all.rds") %>%
  .[["H2000"]] %>% lapply(., function(x) { crs(x) <- "ESRI:31494"; rast(x)}) %>% .[[1]] # for colinearity check

colinvar2000_bv <- correct_colinvar(preds2000_bv, method = c("vif", th = "5"))
env_layer_bv <- terra::wrap(colinvar2000_bv$env_layer)

saveRDS(env_layer_bv, "./1_Inputs/2_Predictors/1_Current/biovars_H2000.rds")

saveRDS(colinvar2000_bv, "./1_Inputs/2_Predictors/1_Current/biovars_H2000_colinvar.rds") # to keep model info

predsrun_bv <- readRDS("./1_Inputs/2_Predictors/1_Current/biovars_all.rds") %>%
  .[[approach]] %>% lapply(., function(x) { crs(x) <- "ESRI:31494"; rast(x)}) %>%
  lapply(., function(x) x %>% select(-colinvar2000_bv$removed_variables) %>% terra::wrap())
saveRDS(predsrun_bv, paste0("./1_Inputs/2_Predictors/1_Current/biovars_", approach, ".rds"))

predsdec_bv <- readRDS("./1_Inputs/2_Predictors/1_Current/biovars_all.rds") %>%
  .[["run10"]] %>% lapply(., function(x) { crs(x) <- "ESRI:31494"; rast(x)}) %>% .[paste0("Y",seq(1975,2015,10))] %>%
  lapply(., function(x) x %>% select(-colinvar2000_bv$removed_variables)  %>% terra::wrap())
saveRDS(predsdec_bv, paste0("./1_Inputs/2_Predictors/1_Current/biovars_dec.rds"))

# Worldclim ====
predictors_wc <- readRDS("./1_Inputs/2_Predictors/1_Current/worldclim_full.rds")[[1]] %>% 
  terra::unwrap(.) %>%
  .[[names(.) %>% str_remove("bio_") %>% as.numeric() %>% sort %>% paste0("bio_", .)]] %>% 
  `names<-`(str_remove(names(.), "_")) %>%
  select(-colinvar2000_bv$removed_variables) %>% 
  terra::wrap()

saveRDS(predictors_wc, "./1_Inputs/2_Predictors/1_Current/worldclim_orig_H2000.rds")

predictors_reproj <- readRDS("./1_Inputs/2_Predictors/1_Current/worldclim_full.rds")[[2]] %>% 
  terra::unwrap(.) %>%
  .[[names(.) %>% str_remove("bio_") %>% as.numeric() %>% sort %>% paste0("bio_", .)]] %>% 
  `names<-`(str_remove(names(.), "_")) %>%
  select(-colinvar2000_bv$removed_variables) %>% 
  terra::wrap()

saveRDS(predictors_reproj, "./1_Inputs/2_Predictors/1_Current/worldclim_reproj_H2000.rds")

# Multivariate indexes =====
preds2000_ind <- readRDS("./1_Inputs/2_Predictors/1_Current/H2000_ind.rds") %>% 
  sapply(., terra::unwrap) %>% 
  .[[1]] %>% `names<-`(str_remove(names(.), "run30_"))

colinvar2000_ind <- correct_colinvar(preds2000_ind %>% 
                                       select(names(.) %>% str_subset(filt, negate = TRUE)), 
                                     method = c("vif", th = "5"))
env_layer_ind <- terra::wrap(colinvar2000_ind$env_layer)

saveRDS(env_layer_ind, "./1_Inputs/2_Predictors/1_Current/ind_H2000.rds")

saveRDS(colinvar2000_ind, "./1_Inputs/2_Predictors/1_Current/ind_H2000_colinvar.rds")

predsrun_ind <- readRDS(paste0("./1_Inputs/2_Predictors/1_Current/", approach, "_ind.rds")) %>%
  sapply(., function(x){
    x <- terra::unwrap(x)
    names(x) <- str_remove(names(x), paste0(approach, "_"))
    x <- x %>% select(names(env_layer_ind %>% terra::unwrap())) %>% terra::wrap()
    return(x)})

saveRDS(predsrun_ind, paste0("./1_Inputs/2_Predictors/1_Current/ind_", approach, ".rds"))

predsdec_ind <- readRDS(paste0("./1_Inputs/2_Predictors/1_Current/", "run10", "_ind.rds")) %>%
  sapply(., function(x){ 
    x <- terra::unwrap(x)
    names(x) <- str_remove(names(x), paste0("run10", "_"))
    x <- x %>% select(names(env_layer_ind %>% terra::unwrap())) %>% terra::wrap()
    return(x)}) %>%
  .[paste0("Y",seq(1975,2015,10))]

saveRDS(predsdec_ind, paste0("./1_Inputs/2_Predictors/1_Current/ind_dec.rds"))

# BOTH ======
preds2000_both <- c(preds2000_ind, preds2000_bv)

colinvar2000_both <- correct_colinvar(preds2000_both %>% 
                                        select(names(.) %>% str_subset(filt, negate = TRUE)), 
                                      method = c("vif", th = "5"))
env_layer_both <- terra::wrap(colinvar2000_both$env_layer)

saveRDS(env_layer_both, "./1_Inputs/2_Predictors/1_Current/both_H2000.rds")

saveRDS(colinvar2000_both, "./1_Inputs/2_Predictors/1_Current/both_H2000_colinvar.rds")

predsrun_ind <- readRDS(paste0("./1_Inputs/2_Predictors/1_Current/", approach, "_ind.rds")) %>%
  sapply(., function(x){ x <- terra::unwrap(x); names(x) <- str_remove(names(x), paste0(approach, "_")); return(x)})

predsrun_bv <- readRDS("./1_Inputs/2_Predictors/1_Current/biovars_all.rds") %>%
  .[[approach]] %>% lapply(., function(x) { crs(x) <- "ESRI:31494"; rast(x)})

predsrun_both <- sapply(names(predsrun_ind), simplify = FALSE, 
                        function(x) c(predsrun_ind[[x]], predsrun_bv[[x]]) %>% 
                          select(names(env_layer_both %>% terra::unwrap())) %>% terra::wrap())

saveRDS(predsrun_both, paste0("./1_Inputs/2_Predictors/1_Current/both_", approach, ".rds"))

predsdec_ind <- readRDS(paste0("./1_Inputs/2_Predictors/1_Current/", approach, "_ind.rds")) %>%
  sapply(., function(x){ x <- terra::unwrap(x); names(x) <- str_remove(names(x), paste0(approach, "_")); return(x)}) %>%
  .[paste0("Y",seq(1975,2015,10))]

predsdec_bv <- readRDS("./1_Inputs/2_Predictors/1_Current/biovars_all.rds") %>%
  .[["run10"]] %>% lapply(., function(x) { crs(x) <- "ESRI:31494"; rast(x)}) %>% .[paste0("Y",seq(1975,2015,10))]

predsdec_both <- sapply(names(predsdec_ind), simplify = FALSE, 
                        function(x) c(predsdec_ind[[x]], predsdec_bv[[x]]) %>% 
                          select(names(env_layer_both %>% terra::unwrap())) %>% terra::wrap())

saveRDS(predsdec_both, paste0("./1_Inputs/2_Predictors/1_Current/both_dec.rds"))
