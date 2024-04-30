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

# ARGUMENTS TO PARSE =======

approach <- c(cmdArg("approach"), cmdArg("a")) # "run10" # options run5, run10, run15, H2000 and worldclim (WC)
preds <- c(cmdArg("p"), cmdArg("preds"), cmdArg("predictors")) # "biovars" # read as BV (bv, biovars) also, and indexes or ind
Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID")) #"centaurea" #
folds <- as.numeric(c(cmdArg("f"), cmdArg("folds")))
wc_type <- c(cmdArg("wc")) #"orig" # reproj
TrimTemp <- ifelse(is.null(c(cmdArg("trimobs"), cmdArg("trim"))), FALSE, TRUE) # trim observations only in during historical period
TrimLon <- ifelse(is.null(c(cmdArg("trimlon"), cmdArg("tlon"), cmdArg("TrimLon"))), FALSE,
                  as.numeric(c(cmdArg("trimlon"), cmdArg("tlon"),cmdArg("TrimLon")))) # trim observations only in during historical period
TrimLat <- ifelse(is.null(c(cmdArg("trimlat"), cmdArg("tlat"), cmdArg("TrimLat"))), FALSE,
                  as.numeric(c(cmdArg("trimlat"), cmdArg("tlat"),cmdArg("TrimLat")))) # trim observations only in during historical period
ThinAlg <- c(cmdArg("thin_alg"), cmdArg("ThinAlg"))
ThinObs <- ifelse(is.null(ThinAlg), FALSE, TRUE)
spat_thr <- as.numeric(c(cmdArg("spt"), cmdArg("spat_thr")))
temp_thr <- as.numeric(c(cmdArg("tt"), cmdArg("temp_thr")))

if(str_detect(preds, "biovars|BV|bv")) {preds <- "biovars"}
if(str_detect(preds, "wc|WC|worldclim")) {preds <- "worldclim"}
if(str_detect(preds, "ind|indexes")) {preds <- "ind"}
if(str_detect(preds, "mix|both")) {preds <- "both"}
if(str_detect(approach, "Hist|hist|2000|h|H")) {approach <- "H2000"}

if(preds == "worldclim") {preds <- paste0(preds, "_", wc_type)}

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

preds_approach <- paste0(preds, "_", approach)

if(approach=="H2000") {
  if(TrimTemp){trim_string <- "_ObsTrim"
  } else {trim_string <- "_NoTrim"}
} else {trim_string <- "_ObsTrim"}

if(is.numeric(TrimLon) & is.numeric(TrimLat)){
  trim_string <- paste0(trim_string, "#TrimLon>",TrimLon, "#TrimLat>",TrimLat)
} else if(is.numeric(TrimLon)){
  trim_string <- paste0(trim_string, "#TrimLon>",TrimLon)
} else if(is.numeric(TrimLat)){
  trim_string <- paste0(trim_string, "#TrimLat>",TrimLat)
}

fit_name <- paste0(preds_approach, trim_string)

if(ThinObs){
  fit_name <- paste0(fit_name, "_ObsThin_",
                     spat_thr/1000, "km",
                     ifelse(length(temp_thr)==0,"", 
                            paste0("_", temp_thr, "Y")))
} else {fit_name <- paste0(fit_name, "_NoThin")}

obs_name <- str_remove(fit_name, paste0(preds, "_"))

print(obs_name)

# DIRECTORIES ==================================================================
Dir.Base <- getwd()

project_directory <- sdm_directory(
  main_dir = Dir.Base,
  projections = c("RCP85", "SSP585"),
  calibration_area = TRUE,
  algorithm = c("all"),
  ensemble = c("meanw"),
  threshold = TRUE,
  return_vector = TRUE
)

# Add folders for the trained models by specie
project_directory %>% str_subset("1_con$") %>%
  str_replace("1_con", paste0("0_models/", Spec4folders)) %>%
  sapply(., function(x) dir.create(x, recursive = TRUE, showWarnings=FALSE)) %>% invisible

dir.create(paste0("./1_Inputs/1_Occurrences/", Spec4folders), 
           recursive = TRUE, showWarnings = FALSE)

# READING DATA =================================================================

## Political Boundaries & DEM ----
CalibArea <- read_sf("./1_Inputs/3_Calibration_area/Kreis_sud_2km.gpkg")
CalibArea_rast <- read_sf("./1_Inputs/3_Calibration_area/vectorized.gpkg")

## Calibration Area ----
CalibArea_terra <- terra::vect(CalibArea)

## Climatic Data ----

preds2use <- readRDS(paste0("./1_Inputs/2_Predictors/1_Current/", preds_approach, ".rds"))
if (is(preds2use, "SpatRaster")) {preds2use <- terra::unwrap(preds2use)} else {preds2use <- lapply(preds2use, terra::unwrap) }

## Occurrences ----

if(file.exists(paste0("./1_Inputs/1_Occurrences/", Spec4folders, "/", obs_name, ".rds"))){
  print("Already thinned...")
  spp_df <- readRDS(paste0("./1_Inputs/1_Occurrences/", Spec4folders, "/", obs_name, ".rds"))[["Train"]]
} else {
  print("Thinning...")
  spp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") # Species data
  
  # MODELLING PROCESS ============================================================
  
  ## Occurrence Subsetting ----
  P_df_full <- spp[spp$Art_wiss == Spec_ID, ] %>%
    mutate(Year = Jahr, 
           x = st_coordinates(geom)[,"X"],
           y = st_coordinates(geom)[,"Y"]) %>% 
    filter(Year <= 2015) %>% 
    arrange(Year, x, y) %>%
    distinct(Year, x, y, .keep_all = TRUE) %>%
    st_intersection(CalibArea_rast) %>%
    st_drop_geometry %>%
    select(x,y,Year) %>% 
    mutate(id = 1:nrow(.)) %>%
    as.data.frame()
  
  P_df_ext <- NULL
  
  if(approach == "H2000"){
    if(TrimTemp){
      P_df <- P_df_full %>%
        filter(Year >= 1970, Year <= 2000)
      P_df_ext <- anti_join(P_df_full,P_df) %>% filter(Year>=1970)
    } else {
      P_df <- P_df_full
    }
  } else if(str_detect(approach, "run")) {
      P_df <- P_df_full %>% filter(Year >= names(preds2use) %>% str_remove("Y") %>% as.numeric %>% min)
  } else if(str_detect(approach, "dec")) {
      P_df <- P_df_full %>% filter(Year >= (names(preds2use) %>% str_remove("Y") %>% as.numeric %>% min ) - 10)
  }
  
  if(is.numeric(TrimLon)){
    P_df <- P_df %>% filter(x<=TrimLon)
  }
  
  if(is.numeric(TrimLat)){
    P_df <- P_df %>% filter(y<=TrimLat)
  }
  
  P_df_val <- P_df
  
  if(ThinObs){
    if (str_detect(approach, "run")) {
      set.seed(42)
      if(ThinAlg == "spt"){
        P_df <- spt_thin(P_df, spat_thr, temp_thr, 20) %>% .[[1]]
      } else {P_df <- spThinAlgMod(P_df, spat_thr, 20) %>% .[[1]]}
    } else if(str_detect(approach, "dec")) {
      P_df <- P_df %>% mutate(Year_i = Year,
                              Year = case_when(Year <= 1975 ~ 1975,
                                               Year > 1975 & Year <= 1985 ~ 1985,
                                               Year > 1985 & Year <= 1995 ~ 1995,
                                               Year > 1995 & Year <= 2005 ~ 2005,
                                               Year > 2005 & Year <= 2015 ~ 2015,
                                               TRUE ~ NA))
      
      P_df_val <- P_df_val %>% mutate(Year_i = Year,
                                      Year = case_when(Year <= 1975 ~ 1975,
                                                       Year > 1975 & Year <= 1985 ~ 1985,
                                                       Year > 1985 & Year <= 1995 ~ 1995,
                                                       Year > 1995 & Year <= 2005 ~ 2005,
                                                       Year > 2005 & Year <= 2015 ~ 2015,
                                                       TRUE ~ NA))
      
      set.seed(42)
      if(ThinAlg == "spt"){
        P_df <- spt_thin(P_df, spat_thr, temp_thr, 20) %>% .[[1]]
      } else {P_df <- spThinAlgMod(P_df, spat_thr, 20) %>% .[[1]]}
    } else if(approach == "H2000"){
      set.seed(42)
      if(ThinAlg == "spt"){
        P_df <- spt_thin(P_df, spat_thr, temp_thr, 20) %>% .[[1]]
      } else {P_df <- spThinAlgMod(P_df, spat_thr, 20) %>% .[[1]]}
    }
  } else if(approach == "dec") {
    P_df <- P_df %>% mutate(Year_i = Year,
                            Year = case_when(Year <= 1975 ~ 1975,
                                             Year > 1975 & Year <= 1985 ~ 1985,
                                             Year > 1985 & Year <= 1995 ~ 1995,
                                             Year > 1995 & Year <= 2005 ~ 2005,
                                             Year > 2005 & Year <= 2015 ~ 2015,
                                             TRUE ~ NA))
    
    P_df_val <- P_df_val %>% mutate(Year_i = Year,
                                    Year = case_when(Year <= 1975 ~ 1975,
                                                     Year > 1975 & Year <= 1985 ~ 1985,
                                                     Year > 1985 & Year <= 1995 ~ 1995,
                                                     Year > 1995 & Year <= 2005 ~ 2005,
                                                     Year > 2005 & Year <= 2015 ~ 2015,
                                                     TRUE ~ NA))
  }
  P_df_val <- anti_join(P_df_val, P_df)

  P_df$Year <- as.character(P_df$Year)
  P_df_val$Year <- as.character(P_df_val$Year)
  P_df$pr_ab <- 1L
  if(nrow(P_df_val)>0) P_df_val$pr_ab <- 1L
  if(!is.null(P_df_ext)) P_df_ext$pr_ab <- 1L
  
  ## RANDOM =====
  set.seed(42)
  sp_part1 <- part_random(
    data = P_df,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = folds)
  )
  data2use <- sp_part1
  
  # Assign decade to the years to match their corresponding env layers
  dec_assign <- function(data, dec) {
    mutate(data, Year = case_when(
      Year <= as.numeric(dec) & Year > (as.numeric(dec) - 10) ~ as.numeric(dec),
      TRUE ~ Year))
  }

  set.seed(42)
  if (str_detect(approach, "run|dec")) {
    years_preds <- names(preds2use) %>% str_remove("Y")
    PseudoA_df <- lapply(years_preds, function(i){
      year_data <- data2use[which(data2use$Year == i),]

      PseudoA_df <- lapply(unique(year_data$.part), function(j){
        pa <- sample_pseudoabs_mod(
          data = data.frame(year_data[which(year_data$.part == j),]),
          x = "x",
          y = "y",
          n = sum(year_data$.part == j),
          method = c('env_const', env = preds2use[[paste0("Y", i)]]),
          rlayer = preds2use[[paste0("Y", i)]],
          calibarea = CalibArea_terra,
        )
        pa$.part <- j
        pa$Year <- i
        return(pa)
      }) %>% bind_rows
    }) %>% bind_rows
  } else {
    PseudoA_df <- lapply(unique(data2use$.part), function(x){
      pa <- sample_pseudoabs_mod(
        data = data.frame(data2use[which(data2use$.part == x),]),
        x = "x",
        y = "y",
        n = sum(data2use$.part == x),
        method = c('env_const', env = preds2use),
        rlayer = preds2use,
        calibarea = CalibArea_terra,
      )
      pa$.part <- x
      return(pa)
    }) %>% bind_rows
  }

  spp_df <- bind_rows(data2use, PseudoA_df)

  if(!is.null(P_df_ext)){
    set.seed(42)
    PseudoA_df_ext <- sample_pseudoabs_mod(
      data = data.frame(P_df_ext),
      x = "x",
      y = "y",
      n = nrow(P_df_ext),
      method = c('env_const', env = preds2use),
      rlayer = preds2use,
      calibarea = CalibArea_terra,
    )
    spp_df_ext <- bind_rows(P_df_ext %>% mutate(Year=as.character(Year)), PseudoA_df_ext)
  } else {spp_df_ext<- NULL}
  
  if(nrow(P_df_val)>0){
    set.seed(42)
    if (str_detect(approach, "run|dec")) {
      years_preds <- names(preds2use) %>% str_remove("Y")
      PseudoA_df_val <- lapply(years_preds, function(i){
        year_data <- P_df_val[which(P_df_val$Year == i),]
        if(nrow(year_data)==0) return()
        pa <- sample_pseudoabs_mod(
          data = data.frame(year_data),
          x = "x",
          y = "y",
          n = nrow(year_data),
          method = c('env_const', env = preds2use[[paste0("Y", i)]]),
          rlayer = preds2use[[paste0("Y", i)]],
          calibarea = CalibArea_terra,
        )
        pa$Year <- i
        return(pa)
      }) %>% bind_rows
    } else {
      PseudoA_df_val <- sample_pseudoabs_mod(
        data = data.frame(P_df_val),
        x = "x",
        y = "y",
        n = nrow(P_df_val),
        method = c('env_const', env = preds2use),
        rlayer = preds2use,
        calibarea = CalibArea_terra,
      )
    }
    
    spp_df_val <- bind_rows(P_df_val %>% mutate(Year=as.character(Year)), PseudoA_df_val)
    saveRDS(list(Train=spp_df, Val=spp_df_val, Ext=spp_df_ext), 
            paste0("./1_Inputs/1_Occurrences/", Spec4folders, "/", obs_name, ".rds"))
  } else {
    saveRDS(list(Train=spp_df, Val=NULL, Ext=spp_df_ext), 
            paste0("./1_Inputs/1_Occurrences/", Spec4folders, "/", obs_name, ".rds"))
  }
}

# EXTRACT ENVIRONMENTAL VALUES =================================================

if(sum(spp_df$pr_ab == 1) >= 10*ifelse(str_detect(approach, "run|dec"),
                                       preds2use[[1]] %>% names %>% length,
                                       preds2use %>% names %>% length)){
  print("Train...")
if (str_detect(approach, "run|dec")) {
  names_preds <- names(preds2use[[1]])
  years_preds <- names(preds2use) %>% str_remove("Y")
  ex_spp <- lapply(years_preds, function(i){
    ex_spp <- sdm_extract(
      data = spp_df %>% filter(Year==i),
      x = "x",
      y = "y",
      env_layer = preds2use[[paste0("Y", i)]], # Raster with environmental variables
      variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
      filter_na = TRUE
    )
  }) %>% bind_rows()
} else {
  names_preds <- names(preds2use)
  ex_spp <- sdm_extract(
    data = spp_df,
    x = "x",
    y = "y",
    env_layer = preds2use, # Raster with environmental variables
    variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
    filter_na = TRUE
  )
}

# FIT AND VALIDATE MODELS ======================================================
## MAX ENTROPY MODELS ========

set.seed(42)
if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/max/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
  max_t <- fit_max(
    data = ex_spp,
    response = "pr_ab",
    predictors = names_preds,
    partition = ".part",
    thr = c("max_sens_spec"),
    clamp = FALSE, # If TRUE, predictors and features are restricted to the range seen during model training. for climate change signal better FALSE?
    classes = "default",
    pred_type = "cloglog",
    regmult = 1
  )
  gc()

  saveRDS(max_t, paste0("./2_Outputs/1_Current/Algorithm/max/0_models/", Spec4folders, "/", fit_name, ".rds"))
} else {
  max_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/max/0_models/", Spec4folders, "/", fit_name, ".rds"))
}

## RANDOM FOREST =========
tune_grid <- expand.grid(mtry = seq(1, 7, 1))

set.seed(42)
if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/raf/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
  rf_t <-tune_raf(
    data = ex_spp,
    response = "pr_ab",
    predictors = names_preds,
    partition = ".part",
    grid = tune_grid,
    thr = "max_sens_spec",
    metric = "SORENSEN",
    n_cores = 1,
  )
  gc()

  saveRDS(rf_t, paste0("./2_Outputs/1_Current/Algorithm/raf/0_models/", Spec4folders, "/", fit_name, ".rds"))
} else {
  rf_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/raf/0_models/", Spec4folders, "/", fit_name, ".rds"))
}

## GAM =====

set.seed(42)
if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/gam/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
  gam_t <- fit_gam(
    data = ex_spp,
    response = "pr_ab",
    predictors = names_preds,
    partition = ".part",
    thr = "max_sens_spec",
  )
  gc()

  saveRDS(gam_t, paste0("./2_Outputs/1_Current/Algorithm/gam/0_models/", Spec4folders, "/", fit_name, ".rds"))
} else {
  gam_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/gam/0_models/", Spec4folders, "/", fit_name, ".rds"))
}

## GAU =====
# Crashing, needs a lot of ram, ignored

# if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/gau/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
#   gau_t <- fit_gau(
#     data = ex_spp,
#     response = "pr_ab",
#     predictors = names_preds,
#     partition = ".part",
#     thr = "max_sens_spec",
#   )
#   saveRDS(gau_t, paste0("./2_Outputs/1_Current/Algorithm/gau/0_models/", Spec4folders, "/", fit_name, ".rds"))
# } else {
#   gau_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/gau/0_models/", Spec4folders, "/", fit_name, ".rds"))
# }
#
# gc()

## GBM =========
tune_grid <-expand.grid(
  n.trees = c(20, 50, 100),
  shrinkage = c(0.1, 0.5, 1),
  n.minobsinnode = c(1, 3, 5, 7, 9)
)

set.seed(42)
if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/gbm/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
  gbm_t <-tune_gbm(
    data = ex_spp,
    response = "pr_ab",
    predictors = names_preds,
    partition = ".part",
    grid = tune_grid,
    thr = "max_sens_spec",
    metric = "SORENSEN",
    n_cores = 1,
  )
  gc()

  saveRDS(gbm_t, paste0("./2_Outputs/1_Current/Algorithm/gbm/0_models/", Spec4folders, "/", fit_name, ".rds"))
} else {
  gbm_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/gbm/0_models/", Spec4folders, "/", fit_name, ".rds"))
}

## GLM =====

set.seed(42)
if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/glm/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
  glm_t <- fit_glm(
    data = ex_spp,
    response = "pr_ab",
    predictors = names_preds,
    partition = ".part",
    thr = "max_sens_spec",
    poly = 2,
    inter_order = 0
  )
  gc()

  saveRDS(glm_t, paste0("./2_Outputs/1_Current/Algorithm/glm/0_models/", Spec4folders, "/", fit_name, ".rds"))
} else {
  glm_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/glm/0_models/", Spec4folders, "/", fit_name, ".rds"))
}

## Neural Net =========
tune_grid <- expand.grid(
  size = c(2, 4, 6, 8, 10),
  decay = c(0.001, 0.05, 0.1, 1, 3, 4, 5, 10)
)

set.seed(42)
if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/net/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
  net_t <-tune_net(
    data = ex_spp,
    response = "pr_ab",
    predictors = names_preds,
    partition = ".part",
    grid = tune_grid,
    thr = "max_sens_spec",
    metric = "SORENSEN",
    n_cores = 1,
  )
  gc()

  saveRDS(net_t, paste0("./2_Outputs/1_Current/Algorithm/net/0_models/", Spec4folders, "/", fit_name, ".rds"))
} else {
  net_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/net/0_models/", Spec4folders, "/", fit_name, ".rds"))
}

## Support Vector Machine =========
tune_grid <- expand.grid(
  C = c(2, 4, 8, 16, 20),
  sigma = c(0.01, 0.1, 0.2, 0.3, 0.4)
)

set.seed(42)
try(
  if(!file.exists(paste0("./2_Outputs/1_Current/Algorithm/svm/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
    svm_t <-tune_svm(
      data = ex_spp,
      response = "pr_ab",
      predictors = names_preds,
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "SORENSEN",
      n_cores = 1,
    )
    gc()

    saveRDS(svm_t, paste0("./2_Outputs/1_Current/Algorithm/svm/0_models/", Spec4folders, "/", fit_name, ".rds"))
  } else {
    svm_t <- readRDS(paste0("./2_Outputs/1_Current/Algorithm/svm/0_models/", Spec4folders, "/", fit_name, ".rds"))
  }
)

# FIT ENSEMBLE ==============

models_ens <- eval(rlang::parse_expr(paste0("list(", ls(pattern = "_t$") %>% 
                                              str_subset("thin", negate = TRUE) %>% 
                                              paste0(., collapse = ", "), ")")))

if(!file.exists(paste0("./2_Outputs/1_Current/Ensemble/meanw/0_models/", Spec4folders, "/", fit_name, ".rds"))) {
  an_ensemble <- fit_ensemble(
    models = models_ens[lengths(models_ens) != 0],
    ens_method = "meanw",
    thr = NULL,
    thr_model = "max_sens_spec",
    metric = "SORENSEN"
  )
  gc()

  saveRDS(an_ensemble, paste0("./2_Outputs/1_Current/Ensemble/meanw/0_models/", Spec4folders, "/", fit_name, ".rds"))
} else {
  an_ensemble <- readRDS(paste0("./2_Outputs/1_Current/Ensemble/meanw/0_models/", Spec4folders, "/", fit_name, ".rds"))
}
} else {
  print("Not enough observations (10x nPred) to train. Exiting...")
}
