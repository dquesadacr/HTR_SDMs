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
  "R.utils",
  "patchwork",
  "ggh4x",
  "fuzzySim",
  "pbapply",
  "brms",
  "lme4",
  "lmerTest",
  "cowplot",
  "tidyr",
  "tidybayes"
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

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

# DIRECTORIES ==================================================================
Dir.Base <- getwd()
setwd(Dir.Base)
Dir.Outputs <- file.path(Dir.Base)
Dir.CasesFolds <- "./"

# READING DATA =================================================================
dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")

MetaEval_ls <- pblapply(uniq_sp, FUN = function(Spec_ID){
  Spec4folders <- str_replace(Spec_ID, " ", "_")

  ## need to be in a specific C..F.. directory
  CF_ls <- lapply(Dir.CasesFolds, FUN = function(Dir.Iter){
      setwd(file.path(Dir.Outputs, Dir.Iter))

    ## Sørensen
    eval_fs <- list.files(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", Spec4folders),
                          pattern = "*.rds$",  full.names = TRUE)
    if(length(eval_fs) == 0){
      eval_all <- NULL
    }else{
      eval_all <- sapply(eval_fs %>% str_subset("eval"),
                         simplify = FALSE, readRDS) %>% bind_rows()
      eval_all$Case <- as.numeric(gsub("\\D", "", unlist(strsplit(Dir.Iter, split = "_"))))[1]
      eval_all$Folds <- as.numeric(gsub("\\D", "", unlist(strsplit(Dir.Iter, split = "_"))))[2]
      eval_all$Species <- Spec_ID

      eval_all <- eval_all %>%
        mutate(Train = Train  %>% str_replace_all(dictio_replace),
               Train = Train %>% str_remove_all("Obs"),
               ObsSubset = ObsSubset %>% str_remove_all("Obs"),
               ObsSubset = ObsSubset %>% str_replace_all(dictio_replace),
               Predictors = Predictors %>% str_replace_all(dictio_replace))  %>%
        separate_wider_delim(Train, delim = "_", too_few = "align_start",
                             names= c("TrainPredictors", "TrainApproach", "TrainTrim", "TrainThin",
                                      "TrainThinDist", "TrainThinTemp")) %>%
        separate_wider_delim(ObsSubset, delim = "_", too_few = "align_start",
                             names= c("EvalObsApproach", "EvalTrim", "EvalThin",
                                      "EvalThinDist", "EvalThinTemp")) %>%
        separate_wider_delim(Predictors, delim = "_", too_few = "align_start",
                             names= c("ExPredictors", "ExPredsApproach")) %>%
        mutate(ExPredsApproach = str_to_sentence(ExPredsApproach)) %>%
        pivot_longer(c(TPR, TNR, SORENSEN, JACCARD, FPB, OR, TSS, AUC, BOYCE, IMAE),
                     names_to = "metric") %>%
        replace_na(list(EvalThinTemp="N/A", EvalThinDist="N/A",
                        TrainThinTemp="N/A", TrainThinDist="N/A")) %>%
        mutate(Posteriori = case_match(Posteriori, TRUE ~ "Post", FALSE ~ "Reg"),
               Posteriori = factor(Posteriori, c("Reg", "Post")),
               TrainApproach = factor(str_to_sentence(TrainApproach),
                                      c("H2k", "Dec", "R10")),
               EvalObsApproach = factor(str_to_sentence(EvalObsApproach),
                                        c("H2k", "Dec", "R10")),
               TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
                                        c("WC", "WC-R", "BV", "MI", "B&I")),
               ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
                                     c("WC", "WC-R", "BV", "MI", "B&I")),
               EvalThinTemp = fct_relevel(EvalThinTemp %>% str_replace_all("Y", ""), c("NA", "0", "1", "5", "10", "100")), #"1Y",
               TrainThinTemp = fct_relevel(TrainThinTemp %>% str_replace_all("Y", ""), c("NA", "0", "1", "5", "10", "100")), #"1Y",
               EvalThinDist = as.numeric(as.character(fct_relevel(EvalThinDist %>% str_replace_all("km", ""), c("NA", "1", "2")))),
               TrainThinDist = as.numeric(as.character(fct_relevel(TrainThinDist %>% str_replace_all("km", ""), c("NA", "1", "2")))),
               Subset = fct_relevel(Subset, c("Train", "Val", "Ext", "Full"))) %>%
        mutate(EvalThinTemp = as.numeric(as.character(fct_relevel(EvalThinTemp %>% str_replace_all("100", "70"), c("NA", "0", "1", "5", "10", "70")))),
               TrainThinTemp = as.numeric(as.character(fct_relevel(TrainThinTemp %>% str_replace_all("100", "70"), c("NA", "0", "1", "5", "10", "70"))))
        )
    }

    ## Schoener's D
    niche_fs <- list.files(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", Spec4folders),
                           pattern = "niche_metrics.rds",  full.names = TRUE)
    if(length(niche_fs) == 0){
      niche_fin <- NULL
    }else{
      niche_metrics <- (niche_fs) %>% readRDS
      dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
      dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")
      niche_fin <- niche_list_process(niche_metrics) %>%
        pivot_longer(cols=c("SchoenerD", "WarrenI", "HellingerDist"))%>%
        filter(str_detect(PredsGRCM, "EC-COS", negate = TRUE), name == "SchoenerD")
      niche_fin <- niche_fin[niche_fin$name == "SchoenerD", ]
      niche_fin$Case <- as.numeric(gsub("\\D", "", unlist(strsplit(Dir.Iter, split = "_"))))[1]
      niche_fin$Folds <- as.numeric(gsub("\\D", "", unlist(strsplit(Dir.Iter, split = "_"))))[2]
      niche_fin$Species <- Spec_ID
    }

    ## return objects
    list(SORENSEN = eval_all,
         SCHOENER = niche_fin)
  })

  ## extract data into return list for each species
  Spec_ls <- list(SORENSEN = NULL,
                  SCHOENER = NULL)
  for(Metr_iter in names(Spec_ls)){
    Iter_ls <- lapply(CF_ls, "[[", Metr_iter)
    Iter_ls <- Iter_ls[!unlist(lapply(Iter_ls, is.null))]
    Iter_df <- do.call(rbind, Iter_ls)
    Spec_ls[[which(names(Spec_ls) == Metr_iter)]] <- Iter_df
  }
  setwd(Dir.Base)
  Spec_ls
})

setwd(Dir.Base)

saveRDS(MetaEval_ls,"MetaEval_ls.rds")

SORENSEN_df <- do.call(rbind, lapply(MetaEval_ls, "[[", "SORENSEN"))

# SORENSEN ANALYSIS ============================================================
## Data Subsetting ----

SORENSEN_df <- SORENSEN_df[SORENSEN_df$metric == "SORENSEN" , ] #& SORENSEN_df$threshold=="max_sens_spec"

SORENSEN_df$EvalThin[SORENSEN_df$EvalThin == "Thin" & is.na(SORENSEN_df$EvalThinTemp)] <- "spThin"
SORENSEN_df$TrainThin[SORENSEN_df$TrainThin == "Thin" & is.na(SORENSEN_df$TrainThinTemp)] <- "spThin"
SORENSEN_model_df <- (SORENSEN_df[, c("value", "thr_value", "ExPredictors",
                                      "EvalObsApproach", "EvalTrim",
                                      "EvalThin", "EvalThinDist", "EvalThinTemp",
                                      "TrainApproach", "TrainTrim",
                                      "TrainThin", "TrainThinDist", "TrainThinTemp",
                                      "Species", "Year",
                                      "Folds", "Case", "n_presences", "Subset", "threshold")])
colnames(SORENSEN_model_df) <- c("SORENSEN", "max_sens_spec", "PREDICTORS",
                                 "Pred_TW", "Pred_TRIM",
                                 "Pred_THIN", "Pred_TD", "Pred_TT",
                                 "Train_TW", "Train_TRIM",
                                 "Train_THIN", "Train_TD", "Train_TT",
                                 "Species", "Year",
                                 "Folds", "Case", "n", "Subset", "threshold")

saveRDS(SORENSEN_model_df, "SORENSEN_model_df.rds")
