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

# READING DATA =================================================================

str2clean <- paste0("_32.nc.rds.xz|_32.rds.xz|evaluation_|ECMWF-|ICTP-|MPI-M-|NCC-|ICHEC-|MOHC-|CLMcom-ETH-|CNRM-|-crCLIM-v1-1|i1p1", "|", 
                    paste0(c("INT","CERFACS-", "-EARTH", "HadGEM2-", "-ESM-LR","ESM1-M"), collapse="|"), "|", paste0(c("DIN63", "MO", "CM4-6"), collapse="|"))

expmts <- c("Hist_RCM", "RCP85")

sapply(expmts, function(expmt){
  print(expmt)
  ifelse(expmt == "Hist_RCM", 
         {expmt_fold <- "1_Current/Hist_RCM"; seq_years <- seq(1965,2005,10); suffix <- "H2000"; run_max <- 30}, 
         {expmt_fold <- paste0("2_Projection/", expmt); seq_years <- seq(2020,2100,10); suffix <- "P20"; run_max <- 20})
  
  ## Climatic Data ----
  set.seed(42)
  # add logic for looping and ordering in rcpXX folders
    
  members <- list.files(paste0("./1_Inputs/2_Predictors/", expmt_fold),
                        pattern = "rds$") %>%
    str_subset("^ind|^bio|^mix|^both|^world", negate = TRUE) %>%
    str_split_i("_bio|_run", 1) %>% unique
  
  sapply(members, FUN = function(y){
    print(y)
    file_name <- y %>% str_remove_all(str2clean) %>% str_split("_", simplify= T) %>% .[,c(1,4,5,3,2,6)] %>% paste0(collapse = "-")

    # Biovars
    preds2000 <- readRDS("./1_Inputs/2_Predictors/1_Current/biovars_H2000.rds")
    
    predsrun <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_biovars.rds")) %>%
      .[[approach]] %>% lapply(., rast)  %>%
      lapply(., function(x) x %>% select(names(preds2000)) %>% terra::wrap())
  
    saveRDS(predsrun, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/biovars_", file_name, "_", approach, ".rds"))
  
    predsdec <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_biovars.rds")) %>%
    .[["run10"]] %>% lapply(., rast) %>% .[paste0("Y",seq_years)] %>% 
      lapply(., function(x) x %>% select(names(preds2000))  %>% terra::wrap())
    saveRDS(predsdec, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/biovars_", file_name,  "_dec.rds"))
  
    predsClim <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_biovars.rds")) %>%
      .[str_subset(names(.), "run10", negate = TRUE)] %>% unlist(recursive = FALSE) %>% lapply(., rast) %>%
      lapply(., function(x) x %>% select(names(preds2000))  %>% terra::wrap())
    names(predsClim) <- names(predsClim) %>% str_split_i("\\.", 2)
    saveRDS(predsClim, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/biovars_", file_name,  "_" , suffix, ".rds"))
  
    # Ind
    preds2000 <- readRDS("./1_Inputs/2_Predictors/1_Current/ind_H2000.rds")
  
    predsrun <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_run10_ind.rds")) %>%
      sapply(., function(x){ x <- terra::unwrap(x); names(x) <- str_remove(names(x), paste0(approach, "_")); x <- x %>% select(names(preds2000)) %>% terra::wrap(); return(x)})
  
    saveRDS(predsrun, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/ind_", file_name, "_", approach, ".rds"))
  
    predsdec <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_run10_ind.rds")) %>%
      sapply(., function(x){ x <- terra::unwrap(x); names(x) <- str_remove(names(x), paste0(approach, "_")); x <- x %>% select(names(preds2000)) %>% terra::wrap(); return(x)}) %>%
      .[paste0("Y",seq_years)]
  
    saveRDS(predsdec, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/ind_", file_name,  "_dec.rds"))
  
    predsClim <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_run", run_max, "_ind.rds")) %>%
      sapply(., function(x){ x <- terra::unwrap(x); names(x) <- str_remove(names(x), paste0(approach, "_")); x <- x %>% select(names(preds2000)) %>% terra::wrap(); return(x)}) #%>% .[paste0("Y",seq(2040,2100,20))]

    saveRDS(predsClim, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/ind_", file_name,  "_", suffix, ".rds"))
    
    # both
    preds2000 <- readRDS("./1_Inputs/2_Predictors/1_Current/both_H2000.rds")

    predsrun_ind <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_run10_ind.rds")) %>%
      sapply(., function(x){ x <- terra::unwrap(x); names(x) <- str_remove(names(x), paste0(approach, "_")); x <- x %>% select(intersect(names(x), names(preds2000))); return(x)})

    predsrun_bio <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_biovars.rds")) %>%
      .[[approach]] %>% lapply(., rast)  %>%
      lapply(., function(x) x %>% select(intersect(names(x), names(preds2000))))

    predsrun <- sapply(names(predsrun_ind), function(x){z <- c(predsrun_ind[[x]], predsrun_bio[[x]]) %>% terra::wrap(); return(z)})

    saveRDS(predsrun, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/both_", file_name, "_", approach, ".rds"))

    predsdec <- sapply(paste0("Y",seq_years), function(x){z <- c(predsrun_ind[[x]], predsrun_bio[[x]]) %>% terra::wrap(); return(z)})

    saveRDS(predsdec, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/both_", file_name,  "_dec.rds"))

    predsClim_bio <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_biovars.rds")) %>%
      .[str_subset(names(.), "run10", negate = TRUE)] %>% unlist(recursive = FALSE) %>% lapply(., rast) %>%
      lapply(., function(x) x %>% select(intersect(names(x), names(preds2000))))
    names(predsClim_bio) <- names(predsClim_bio) %>% str_split_i("\\.", 2)

    predsClim_ind <- readRDS(paste0("./1_Inputs/2_Predictors/", expmt_fold, "/", y, "_run", run_max, "_ind.rds")) %>%
      sapply(., function(x){ x <- terra::unwrap(x); names(x) <- str_remove(names(x), paste0(approach, "_")); x <- x %>% select(intersect(names(x), names(preds2000))); return(x)})

    predsClim <- sapply(names(predsClim_ind), function(x){z <- c(predsClim_ind[[x]], predsClim_bio[[x]]) %>% terra::wrap(); return(z)})
    saveRDS(predsClim, paste0("./1_Inputs/2_Predictors/", expmt_fold, "/both_", file_name,  "_", suffix, ".rds"))
  })
})
