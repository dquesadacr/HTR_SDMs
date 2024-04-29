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
preds <- c(cmdArg("p"), cmdArg("preds"), cmdArg("predictors"))
plot_bool <- ifelse(is.null(c(cmdArg("plot"), cmdArg("pl"))), FALSE, TRUE)
Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID"))
ow_bool <- ifelse(is.null(c(cmdArg("ow"), cmdArg("o"), cmdArg("overwrite"))), FALSE, TRUE)

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

ens_list <- sapply(list.files(paste0("./2_Outputs/1_Current/Ensemble/meanw/0_models/", Spec4folders),
                              paste0("^", preds,"*"), full.names = TRUE), simplify = FALSE, readRDS)
names(ens_list) <- names(ens_list) %>% str_split_i("/", -1) %>% 
  str_remove(".rds")

preds_list <- sapply(list.files(paste0("./1_Inputs/2_Predictors/1_Current"),
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

preds_full <- str_remove_all(names(preds_list),
                             (paste0("_", str_split_i(names(preds_list), "_", -1) %>% 
                                       unique))) %>% unique %>% paste0(., "_")
if(length(preds_full)>1) {preds_full <- paste0(preds_full, collapse = "|")}

TV_list <- sapply(list.files(paste0("./1_Inputs/1_Occurrences/", Spec4folders),
                             ".rds$", full.names = TRUE), simplify = FALSE, readRDS)

names(TV_list) <- names(TV_list) %>% str_split_i("/", -1) %>% str_remove_all(".rds")

# PREDICT ====

plot_preds <- function(ind_p, ens_train, preds, period, posteriori = FALSE){
  if(posteriori) {
    ens_type <- "predict_posteriori_"
    data2plot <- ind_p
  } else {
    ens_type <- "predict_"
    data2plot <- ind_p$meanw
  }
  dir.create(paste0("./2_Outputs/1_Current/Ensemble/meanw/3_plots/",
                    Spec4folders, "/", ens_type,
                    ens_train, "/", preds),
             recursive = TRUE, showWarnings = FALSE)

  plot_ind <- ggplot() +
  geom_spatraster(data = data2plot$max_sens_spec) +
  scale_fill_viridis_d(direction = 1, name ="", begin = 0.2, end = 0.8, na.value = NA, breaks=c("TRUE", "FALSE")) +
  new_scale_fill() +
  geom_spatraster(data = data2plot$meanw) +
  scale_fill_whitebox_c(palette = "soft",na.value = NA,
                        name="", direction=-1, limits=c(0,1)) +
  facet_wrap(~lyr) +
  theme_light(base_size = 10) +
  labs(title= period %>% str_remove("Y")) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color= "black"),
        legend.position = "bottom",
        legend.key.width = unit(7.5, "mm"),
        legend.key.height = unit(3, "mm"),
        plot.title = element_text(hjust = 0.5))

  ggsave(plot_ind, filename=paste0("./2_Outputs/1_Current/Ensemble/meanw/3_plots/",
                                   Spec4folders, "/", ens_type, ens_train, "/",
                                   preds, "/", period %>% str_remove("Y"), ".png"),
         width=15, height=8, units="cm")
}

gif_preds <- function(ens_train, preds, posteriori = FALSE){
  if(posteriori) {ens_type <- "predict_posteriori_"} else {ens_type <- "predict_"}

  if(!str_detect(preds,"H2000")){
    system(paste0("convert -delay 20 -loop 0 ./2_Outputs/1_Current/Ensemble/meanw/3_plots/",
                  Spec4folders, "/", ens_type, ens_train, "/", preds, "/*png ",
                  "./2_Outputs/1_Current/Ensemble/meanw/3_plots/",
                  Spec4folders, "/", ens_type, ens_train, "/", preds, "/anim.gif"))
  }
}

ens_predict <- function(ensemble_list, predictors_list, plot = TRUE, overwrite = FALSE) {
  if(!overwrite & file.exists(paste0("./2_Outputs/1_Current/Ensemble/meanw/1_con/", 
                        Spec4folders, "/preds_ens_", preds, ".rds"))){
    # if(!overwrite) {
      pred_all <- readRDS(paste0("./2_Outputs/1_Current/Ensemble/meanw/1_con/",
                                 Spec4folders, "/preds_ens_", preds, ".rds"))
    # }
  } else {
    pred_all <- sapply(names(ensemble_list), simplify = FALSE, function(x){
      per_ens <- sapply(names(predictors_list), simplify = FALSE, function(y){
        per_period <- sapply(names(predictors_list[[y]]), simplify = FALSE, function(z){
          print(paste0(z, " -- ", y, " -- ", x))
          pred <- sdm_predict(
            models = ensemble_list[[x]],
            pred = predictors_list[[y]][[z]],
            thr = "max_sens_spec",
            con_thr = FALSE,
            predict_area = NULL
          )
          if(plot) {plot_preds(ind_p = pred, ens_train = x, preds = y, period = z)}
          pred <- sapply(pred, simplify = FALSE, function(x) {x <- terra::wrap(x); return(x)})
          return(pred)
        })
        if(plot) {gif_preds(ens_train = x, preds = y)}
        return(per_period)
      })
    })
    # Saving the full list with con and bin in con
    dir.create(paste0("./2_Outputs/1_Current/Ensemble/meanw/1_con/", Spec4folders),
              recursive = TRUE, showWarnings = FALSE)
    saveRDS(pred_all, paste0("./2_Outputs/1_Current/Ensemble/meanw/1_con/", Spec4folders, "/preds_ens_", preds, ".rds"))
  }
  return(pred_all)
}

unwrap_lists <- function(x){
  if (is.list(x)) {x <- sapply(x, unwrap_lists, simplify = FALSE)
  } else {x <- terra::unwrap(x)}
}

pred_all <- ens_predict(ens_list, preds_list, plot = plot_bool, overwrite = FALSE) %>% unwrap_lists()

# ENSENSEMBLE EVALUATION ========

eval_df <- function(predictions_list, train_val_list, posteriori = FALSE, overwrite = ow_bool){
  if(posteriori) {ens_type <- "eval_posteriori_"} else {ens_type <- "eval_regular_"}
  dir.create(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", Spec4folders),
             recursive = TRUE, showWarnings = FALSE)
  if(!overwrite & file.exists(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                        Spec4folders, "/", ens_type, preds, ".rds"))){
    eval_all <- readRDS(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", 
                               Spec4folders, "/", ens_type, preds, ".rds"))
  } else {
  eval_all <- sapply(names(predictions_list), simplify = FALSE, function(x){
    per_trama <- sapply(names(predictions_list[[x]]), simplify = FALSE, function(w){

    Trama <- FALSE
      if(str_remove_all(w, preds_full) == (str_remove_all(x, preds_full) %>% str_split_i("_", 1)) ){
        tv_names <- str_subset(names(train_val_list),
                                   str_remove_all(x, ifelse(str_detect(preds_full, "|"), preds_full, w)))
      } else {
        tv_names<- names(train_val_list) %>% str_subset("NoThin") %>%
          str_subset("H2000_Obs", negate=TRUE) %>% str_subset(str_remove_all(w, preds_full) )
        Trama <- TRUE
      }

      per_ens <- sapply(tv_names,simplify = FALSE, function(y){
        per_subset <- sapply(names(train_val_list[[y]]), simplify = FALSE, function(v){
          per_period <- sapply(names(predictions_list[[x]][[w]]), simplify = FALSE, function(z){
          if(z == "H2000"){
            ex_spp <- train_val_list[[y]][[v]]
            if(is.null(ex_spp)) {return()}
          } else {
            if(is.null(train_val_list[[y]][[v]])) {return()}
            ex_spp <- train_val_list[[y]][[v]] %>% filter(Year == z %>% str_remove("Y"))
            if(nrow(ex_spp) == 0 | nrow(ex_spp %>% filter(pr_ab == 0)) == 0) {return()}
          }
          predict_suit <- sdm_extract(
            data = ex_spp,
            x = "x",
            y = "y",
            env_layer = predictions_list[[x]][[w]][[z]]$meanw$meanw,
            variables = NULL,
            filter_na = TRUE
          )
          presence_num <- predict_suit %>% filter(pr_ab == 1) %>% select(meanw) %>% .[[1]] %>% as.numeric
          absence_num <- predict_suit %>% filter(pr_ab == 0) %>% select(meanw) %>% .[[1]] %>% as.numeric

          eval_ens <- sdm_eval(
            p = presence_num,
            a = absence_num
          )

          eval_ens$Year <- z %>% str_remove("Y")
          print(c(x,y,z,w,v))
          return(eval_ens)
        }) %>% bind_rows()
        per_period$Predictors <- w
        per_period$ObsSubset <- y
        per_period$Subset <- v
        if(Trama) {per_period$Subset <- "Full"}
        per_period$Train <- x
        per_period$Posteriori <- posteriori
        return(per_period)
      }) %>% bind_rows()
    }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
  saveRDS(eval_all, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", Spec4folders, "/", ens_type, preds, ".rds"))
  }
  return(eval_all)
}

eval_all <- eval_df(pred_all, TV_list, overwrite=ow_bool)

