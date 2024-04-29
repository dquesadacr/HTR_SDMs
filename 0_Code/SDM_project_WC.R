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
preds <- "wc"
plot_bool <- ifelse(is.null(c(cmdArg("plot"), cmdArg("pl"))), FALSE, TRUE)
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

ens_list <- sapply(list.files(paste0("./2_Outputs/1_Current/Ensemble/meanw/0_models/", Spec4folders),
                              paste0("^", preds,"*"), full.names = TRUE), simplify = FALSE, readRDS)
names(ens_list) <- names(ens_list) %>% str_split_i("/", -1) %>% str_remove(".rds")

#expmt <- "RCP85"
expmt <- "SSP585" # experiment for WC

preds_list <- sapply(list.files(paste0("./1_Inputs/2_Predictors/2_Projection/", expmt),
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

plot_preds <- function(ind_p, ens_train, preds, period, posteriori = FALSE){
  if(posteriori) {
    ens_type <- "predict_posteriori_"
    data2plot <- ind_p
  } else {
    ens_type <- "predict_"
    data2plot <- ind_p$meanw
  }
  dir.create(paste0("./2_Outputs/2_Projection/", expmt, "/Ensemble/meanw/3_plots/",
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

  ggsave(plot_ind, filename=paste0("./2_Outputs/2_Projection/", expmt, "/Ensemble/meanw/3_plots/",
                                   Spec4folders, "/", ens_type, ens_train, "/",
                                   preds, "/", period %>% str_remove("Y"), ".png"),
         width=15, height=8, units="cm")
}

gif_preds <- function(ens_train, preds, posteriori = FALSE){
  if(posteriori) {ens_type <- "predict_posteriori_"} else {ens_type <- "predict_"}

  if(!str_detect(preds,"H2000")){
    system(paste0("convert -delay 20 -loop 0 ./2_Outputs/2_Projection/", expmt, "/Ensemble/meanw/3_plots/",
                  Spec4folders, "/", ens_type, ens_train, "/", preds, "/*png ",
                  "./2_Outputs/2_Projection/", expmt, "/Ensemble/meanw/3_plots/",
                  Spec4folders, "/", ens_type, ens_train, "/", preds, "/anim.gif"))
  }
}

ens_predict <- function(ensemble_list, predictors_list, plot = TRUE) {
  pred_all <- sapply(names(ensemble_list), simplify = FALSE, function(x){
    per_ens <- sapply(names(predictors_list), simplify = FALSE, function(y){
      per_period <- sapply(names(predictors_list[[y]]), simplify = FALSE, function(z){
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
  dir.create(paste0("./2_Outputs/2_Projection/", expmt, "/Ensemble/meanw/1_con/", Spec4folders),
             recursive = TRUE, showWarnings = FALSE)
  saveRDS(pred_all, paste0("./2_Outputs/2_Projection/", expmt, "/Ensemble/meanw/1_con/", Spec4folders, "/preds_ens_", preds, ".rds"))
  return(pred_all)
}

unwrap_lists <- function(x){
  if (is.list(x)) {x <- sapply(x, unwrap_lists, simplify = FALSE)
  } else {x <- terra::unwrap(x)}
}

pred_all <- ens_predict(ens_list, preds_list, plot = plot_bool) %>% unwrap_lists()
