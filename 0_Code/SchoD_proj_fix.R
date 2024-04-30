
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
  "R.utils",
  "patchwork",
  "ggh4x",
  "fuzzySim"
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

Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID"))

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")

pred_all <- sapply(list.files(paste0("./2_Outputs/1_Current/Ensemble/meanw/1_con/", Spec4folders),
                              pattern = "*.rds$",  full.names = TRUE),
                   simplify = FALSE, readRDS) %>% unwrap_lists()

names(pred_all) <- names(pred_all) %>% str_split_i("/", -1) %>% str_remove_all(".rds|preds_")

pred_hist <- sapply(list.files(paste0("./2_Outputs/1_Current/Hist_RCM/Ensemble/meanw/1_con/", Spec4folders),
                               pattern = "*.rds$",  full.names = TRUE) %>% str_subset("ind.rds|both.rds|biovars.rds"),
                    simplify = FALSE, readRDS) %>% unwrap_lists()

names(pred_hist) <- names(pred_hist) %>% str_split_i("/", -1) %>% str_remove_all(".rds|preds_")


if(!file.exists(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                      Spec4folders, "/niche_metrics.rds"))){


mss_predictions <- binPixels2vec(pred_all) # originally with binary maps, but testing with continous
mss_hist <- binPixels2vec(pred_hist)

niche_metrics <- sapply(names(mss_hist), simplify = FALSE, function(z){
  sapply(names(mss_hist[[z]]), simplify = FALSE, function(y){
    sapply(names(mss_hist[[z]][[y]]), simplify = FALSE, function(x){
      sapply(names(mss_hist[[z]][[y]][[x]]), simplify = FALSE, function(w){
        sapply(names(mss_hist[[z]][[y]][[x]][[w]]), simplify = FALSE, function(v){
          sapply(names(mss_hist[[z]][[y]][[x]][[w]][[v]]), simplify = FALSE, function(u){
            z2 <- paste0("ens_", z %>% str_split_i("_", 2))
            ifelse(str_split_i(x,"_", 3) == "H2000", w2 <- "H2000", w2 <- w)
            # print(paste(z,z2,y,x,w, w2, v,u, sep = " - "))
            p1 <- mss_predictions[[z2]][[y]][[paste0(x %>% str_split_i("_", 1), "_", x %>% str_split_i("_", 3))]][[w2]][[v]][[u]]
            p2 <- mss_hist[[z]][[y]][[x]][[w]][[v]][[u]]
            
            if(!is.null(p1)){modOverlap(p1,p2)} # & !is.null(p2)
          })
        })
      })
    })
  })
})
gc()

saveRDS(mss_predictions, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                                Spec4folders, "/mss_predictions.rds"))

saveRDS(mss_hist, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                         Spec4folders, "/mss_hist.rds"))

saveRDS(niche_metrics, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                              Spec4folders, "/niche_metrics.rds"))
}

if(!file.exists(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                       Spec4folders, "/niche_metrics_cont.rds"))){
mss_predictions_cont <- contPixels2vec(pred_all)
mss_hist_cont <- contPixels2vec(pred_hist)

niche_metrics_cont <- sapply(names(mss_hist_cont), simplify = FALSE, function(z){
  sapply(names(mss_hist_cont[[z]]), simplify = FALSE, function(y){
    sapply(names(mss_hist_cont[[z]][[y]]), simplify = FALSE, function(x){
      sapply(names(mss_hist_cont[[z]][[y]][[x]]), simplify = FALSE, function(w){
        sapply(names(mss_hist_cont[[z]][[y]][[x]][[w]]), simplify = FALSE, function(v){
          sapply(names(mss_hist_cont[[z]][[y]][[x]][[w]][[v]]), simplify = FALSE, function(u){
            z2 <- paste0("ens_", z %>% str_split_i("_", 2))
            ifelse(str_split_i(x,"_", 3) == "H2000", w2 <- "H2000", w2 <- w)
            p1 <- mss_predictions_cont[[z2]][[y]][[paste0(x %>% str_split_i("_", 1), "_", x %>% str_split_i("_", 3))]][[w2]][[v]][[u]]
            p2 <- mss_hist_cont[[z]][[y]][[x]][[w]][[v]][[u]]
            
            if(!is.null(p1)){modOverlap(p1,p2)} # & !is.null(p2)
          })
        })
      })
    })
  })
})
gc()

niche_fin_cont <- niche_list_process(niche_metrics_cont) %>%
  pivot_longer(cols=c("SchoenerD", "WarrenI", "HellingerDist")) %>%
  filter(str_detect(PredsGRCM, "EC-COS", negate = TRUE), name == "SchoenerD")

saveRDS(mss_predictions_cont, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                                     Spec4folders, "/mss_predictions_cont.rds"))

saveRDS(mss_hist_cont, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                              Spec4folders, "/mss_hist_cont.rds"))

saveRDS(niche_metrics_cont, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                                   Spec4folders, "/niche_metrics_cont.rds"))
}

dictio_replace_reverse <- c("N/A"=NA, "70Y"="100Y", "H2k"="H2000", " km"="km",
                            "R10"="run10", "Dec"="dec", "WC$"="worldclim_orig",
                            "WC-R"="worldclim_reproj", "MI"="ind", "B&I"="both", 
                            "BV"="biovars", "^Thin"="ObsThin", "^Trim"="ObsTrim")

# Proj
if(!file.exists(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                       Spec4folders, "/perc_elev.rds"))){

pred_rcp <- sapply(c(list.files(paste0("./2_Outputs/2_Projection/RCP85/Ensemble/meanw/1_con/", Spec4folders),
                                pattern = "*.rds$",  full.names = TRUE) %>% str_subset("ind.rds|both.rds|biovars.rds"),
                     list.files(paste0("./2_Outputs/2_Projection/SSP585/Ensemble/meanw/1_con/", Spec4folders),
                                pattern = "*.rds$",  full.names = TRUE) %>% str_subset("ind.rds|both.rds|biovars.rds|clim.rds")),
                   simplify = FALSE, readRDS) %>% unwrap_lists()

names(pred_rcp) <- names(pred_rcp) %>% str_split_i("/", -1) %>% str_remove_all(".rds|preds_")

pixels_total <- as.vector(pred_rcp$ens_biovars$biovars_dec_ObsTrim_NoThin$`biovars_MPI-COS-v1-r1-rcp85-det_dec`$Y2020$meanw$max_sens_spec) %>%
  na.omit() %>% length()

DEM <- rast("./1_Inputs/3_Calibration_area/topo.tif") # Digital Elevation Model
crs(DEM) <- "ESRI:31494"
DEM <- crop(DEM, ext(pred_rcp$ens_biovars$biovars_dec_ObsTrim_NoThin$`biovars_MPI-COS-v1-r1-rcp85-det_dec`$Y2020$meanw$max_sens_spec))

mss_rcp_topoPix <- binPixels2topoPix(pred_rcp)
mss_pred_topoPix <- binPixels2topoPix(pred_all)
mss_hist_topoPix <- binPixels2topoPix(pred_hist)

mss_test3 <- mss_hist_topoPix %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  discard(is.null)%>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(Pred = rownames(.)) %>%
  separate_wider_delim(cols = Pred, delim = "+",
                       names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric", "Var"),
                       too_few = "align_start") %>%
  as_tibble() %>%
  pivot_wider(names_from = Var, values_from = V1) %>%
  mutate(PredictApproach = str_replace_all(PredictApproach, c("H2000"="P20")))

mss_test2 <- mss_pred_topoPix %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  discard(is.null)%>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(Pred = rownames(.)) %>%
  separate_wider_delim(cols = Pred, delim = "+",
                       names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric", "Var"),
                       too_few = "align_start") %>%
  as_tibble() %>%
  pivot_wider(names_from = Var, values_from = V1) %>%
  mutate(PredictApproach = str_replace_all(PredictApproach, c("_"="_Train_","H2000"="P20")))

mss_test <- mss_rcp_topoPix %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  list_flatten(name_spec = "{outer}+{inner}") %>%
  discard(is.null) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(Pred = rownames(.)) %>%
  separate_wider_delim(cols = Pred, delim = "+", 
                       names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric", "Var"), 
                       too_few = "align_start") %>% 
  as_tibble() %>%
  pivot_wider(names_from = Var, values_from = V1) %>%
  rbind(mss_test2, mss_test3, .) %>%
  mutate(Elev = ifelse(is.nan(Elev),NA,Elev),
         Pixels = Pixsum / pixels_total,
         ModelOrig = TrainPredictors,
         Model = TrainPredictors %>% str_replace_all(c(dictio_replace, "Obs"="", dictio_replace_preds)),
         TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
         TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
         PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
  drop_na(Elev) %>%
  separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach")) %>%
  separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
  separate(TrainApproach, sep = "-", c("TrainApproachPeriod", "TrainApproachObsTrim")) %>%
  separate(PredictApproach, sep = "_", 
           c("ExPredictors", "PredsGRCM", "PredsPeriod")) %>%
  separate(PredsGRCM, sep = "-", 
           c("PredsGCM", "PredsRCM", "PredsRCMver", "PredsRealisation", "PredsExp", "PredsDetSto")) %>%
  mutate(TrainApproachObsTrim = case_when(
    TrainApproachPeriod == "H2k" ~ tidyr::replace_na(TrainApproachObsTrim, "NoTrim"),
    TRUE ~ tidyr::replace_na(TrainApproachObsTrim, "ObsTrim")),
    TrainApproachObsTrim = paste0("Train",TrainApproachObsTrim),
    TrainPredsType = case_when(
      !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
      TRUE ~ NA_character_),
    PredsGRCM = paste0(PredsGCM, "-", PredsRCM, "-", PredsRCMver, "-", PredsRealisation),
    PredsExp = case_when(str_detect(PredsExp, "hist") ~ str_to_sentence(PredsExp) %>% str_remove_all("orical"),
                         TRUE ~ str_to_upper(PredsExp)),
    TrainApproachPeriod = factor(
      paste0("Train", str_to_sentence(TrainApproachPeriod)),
      paste0("Train", c("H2k", "Dec", "R10"))),
    PredsPeriod = str_replace_all(PredsPeriod, "rcp85", "P20"),
    PredsPeriod = factor(paste0("Pred", str_to_sentence(PredsPeriod)),
                         paste0("Pred", c("H2k", "P20", "Dec", "R10"))),
    TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
                             c("WC", "BV", "MI", "B&I")),
    ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
                          c("WC", "BV", "MI", "B&I")),
    Year = str_split_i(Year,"\\.", 1) %>% str_remove_all("Y|H") %>% as.numeric()) %>%
  mutate(PredsPeriod=replace_na(PredsPeriod, "PredP20"))

saveRDS(mss_test, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                         Spec4folders, "/perc_elev.rds"))
}
