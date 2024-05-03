
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

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique %>% sort %>% str_subset("rchis", negate=TRUE)

Spec4folders_all <- str_replace(uniq_sp, " ", "_")

dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")

cont_niche_df <- list.files(pattern = "niche_metrics_cont.rds", recursive = TRUE) %>%
  str_subset("mascula", negate = TRUE) %>%
  sapply(., simplify = FALSE, FUN=function(x){
    niche_metrics <- readRDS(x)
    niche_fin <- niche_list_process(niche_metrics) %>%
      pivot_longer(cols=c("SchoenerD", "WarrenI", "HellingerDist")) %>%
      filter(str_detect(PredsGRCM, "EC-COS", negate = TRUE), name == "SchoenerD")
    niche_fin$Species <- str_split_i(x, "/", 5) %>% str_replace_all("_", "\n")
    return(niche_fin)
  }) %>% bind_rows()

niche_filt <- cont_niche_df %>% ungroup %>%
  filter(TrainPredictors %in% c("B&I"),
         TrainApproach != "Dec" & PredsPeriod != "Dec") %>%
  mutate(TrainThinDist = fct_relabel(TrainThinDist, ~ gsub("N/A", "0 km", .x)),
         TrainThinTemp = if_else(TrainThinDist == "0 km",
                                 if_else(TrainThinTemp=="N/A",
                                         fct_relabel(TrainThinTemp, ~ gsub("N/A", "", .x)),
                                         TrainThinTemp),
                                 fct_relabel(TrainThinTemp, ~ gsub("N/A", "spThin", .x))))

h2k_df <- niche_filt %>%
  filter(TrainApproach == "H2k",
         TrainThinTemp %in% c("","spThin"),
         TrainTrim != "Trim")

r10_df <- niche_filt %>%
  filter(TrainApproach == "R10",
         TrainThinTemp == "0Y" & TrainThinDist == "1 km" |
         TrainThinTemp == "5Y" & TrainThinDist == "1 km" |
         TrainThinTemp == "0Y" & TrainThinDist == "2 km" |
         TrainThinTemp == "5Y" & TrainThinDist == "2 km" )

niche_fin <- bind_rows(h2k_df, r10_df)

niche_plot1 <- ggplot(niche_fin,
                      aes(x=Species,
                          y=value, color=Species)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_manual(values = c("#7d5690", "#34acf8", "#6acc68", "#f7b036", "#ef9fef",
                                "#758796", "#1929C8", "#C42503", "#3f822e", "#C1C1C1", "#8F0000"),
                     name="") +
  facet_nested("Predict"+PredsPeriod~TrainApproach+TrainThinDist+TrainThinTemp,
               scales = "free_x", space = "free_x",
               nest_line = element_line(color="black", linewidth = .25)) +
  theme_light(base_size = 10) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", margin = margin(0.5,0.5,0,0, unit = "mm")),
    legend.position = "bottom",
    legend.box.margin = margin(-3, 0, 0.5, 0, unit = "mm"),
    panel.spacing = unit(1, "mm"),
    plot.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.text = element_text(margin=margin(0, 1, 0, 0, unit = "mm"), face = "italic"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing.x = unit(2, "mm"),
    legend.spacing.y = unit(0, "mm"),
  )

####

pred_fin <- sapply(Spec4folders_all, simplify=FALSE, FUN=function(Spec4folders){
  print(Spec4folders)
  pred_all <- sapply(list.files(paste0("./2_Outputs/1_Current/Ensemble/meanw/1_con/", Spec4folders),
                                pattern = "*.rds$",  full.names = TRUE),
                    simplify = FALSE, readRDS) %>% unwrap_lists()

  names(pred_all) <- names(pred_all) %>% str_split_i("/", -1) %>% str_remove_all(".rds|preds_")

  pred_all_latlon <- sapply(list.files(paste0("../../LatLon/C1_F10/2_Outputs/1_Current/Ensemble/meanw/1_con/", Spec4folders),
                                pattern = "*.rds$",  full.names = TRUE),
                    simplify = FALSE, readRDS) %>% unwrap_lists()

  names(pred_all_latlon) <- names(pred_all_latlon) %>% str_split_i("/", -1) %>% str_remove_all(".rds|preds_")

  mss_predictions <- contPixels2vec(pred_all) # originally with binary maps, but testing with continous
  mss_predictions_latlon <- contPixels2vec(pred_all_latlon)

  niche_metrics <- sapply(names(mss_predictions_latlon), simplify = FALSE, function(z){
    sapply(names(mss_predictions_latlon[[z]]), simplify = FALSE, function(y){
      sapply(names(mss_predictions_latlon[[z]][[y]]), simplify = FALSE, function(x){
        sapply(names(mss_predictions_latlon[[z]][[y]][[x]]), simplify = FALSE, function(w){
          sapply(names(mss_predictions_latlon[[z]][[y]][[x]][[w]]), simplify = FALSE, function(v){
            sapply(names(mss_predictions_latlon[[z]][[y]][[x]][[w]][[v]]), simplify = FALSE, function(u){
              z2 <- paste0("ens_", z %>% str_split_i("_", 2))
              p1 <- mss_predictions[[z2]][[y %>% str_remove_all("#TrimLon>4615000#TrimLat>5625000|#TrimLon>4615000|#TrimLat>5625000")]][[x]][[w]][[v]][[u]]
              p2 <- mss_predictions_latlon[[z]][[y]][[x]][[w]][[v]][[u]]

              if(!is.null(p1)){modOverlap(p1,p2)} # & !is.null(p2)
            })
          })
        })
      })
    })
  })
  gc()
  return(niche_metrics)
})

saveRDS(pred_fin, paste0("./sp_trans.rds"))

# pred_fin <- readRDS("./sp_trans.rds")

pred_df <- sapply(names(pred_fin) %>% str_subset("Dianthu", negate = TRUE), 
                  simplify = FALSE, FUN= function(spec){
  niche_fin_cont <- niche_list_process(pred_fin[[spec]]) %>%
    pivot_longer(cols=c("SchoenerD", "WarrenI", "HellingerDist")) %>%
    filter(str_detect(PredsGRCM, "EC-COS", negate = TRUE), name == "SchoenerD")
  niche_fin_cont$Species <- spec %>% str_replace("_", " ")
  niche_fin_cont <- niche_fin_cont %>% filter(PredsGCM==TrainApproach)
  return(niche_fin_cont)
}) %>% bind_rows

pred_df <- pred_df %>%
  mutate(TrainPredictors = str_replace_all(TrainPredictors, "WC", "WC-R"),
         TrainThinDist = fct_relabel(TrainThinDist, ~ gsub("N/A", "0 km", .x)),
         TrainThinTemp = if_else(TrainThinDist == "0 km", 
                                 if_else(TrainThinTemp=="N/A", 
                                         fct_relabel(TrainThinTemp, ~ gsub("N/A", "", .x)),
                                         TrainThinTemp),
                                 fct_relabel(TrainThinTemp, ~ gsub("N/A", "spThin", .x))))

pred_df$Species <- factor(pred_df$Species %>% str_replace_all(" ", "\n"), 
                          uniq_sp %>% str_replace_all(" ", "\n"))

niche_plot2 <- ggplot(pred_df,
                      aes(x=Species,
                          y=value, color=Species)) +
  geom_boxplot(outlier.size = 0.25, size=0.4) +
  scale_color_manual(values = c("#7d5690", "#34acf8", "#6acc68", "#f7b036", "#ef9fef",
                                "#758796", "#1929C8", "#C42503", "#3f822e", "#C1C1C1", "#8F0000"),
                     breaks = uniq_sp %>% str_replace_all(" ", "\n"),
                     name="") +
  scale_x_discrete(drop=FALSE) +
  facet_nested("  "+"1 km"~TrainApproach+TrainPredictors+TrainThinTemp,
               nest_line = element_line(color="black", linewidth = .25)) +
  theme_light(base_size = 10) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", margin = margin(0.4,0,0.4,0, unit = "mm")),
    legend.position = "bottom",
    legend.box.margin = margin(-3, 0, 0.5, 0, unit = "mm"),
    panel.spacing = unit(1, "mm"),
    plot.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.text = element_text(margin=margin(0, 1, 0, 0, unit = "mm"), face = "italic"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing.x = unit(2, "mm"),
    legend.spacing.y = unit(0, "mm"),
  )


# ggsave(plot= niche_plot2, filename=paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", "spat_trans.pdf"),
#        width=180, height=60, units="mm")

####

p1_tot <- ggarrange(plotlist = list(
  niche_plot1 + 
    facet_nested("Predict"+PredsPeriod~TrainApproach+TrainThinDist+TrainThinTemp,
                 nest_line = element_line(color="black", linewidth = .25)) +
    scale_y_continuous(breaks = c(0.7,0.8,0.9)) +
    theme(legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
          plot.margin = margin(0, 0, 0.5, 0, unit = "mm")),
  niche_plot2 + theme(legend.box.margin = margin(0, 0, 0, 0, unit = "mm"))),
                    labels =c("(a)", "(b)"),
                    label.x = -0.01,common.legend = TRUE,
                    label.y = 1.01,
                    legend = "bottom",
                    nrow = 2,
                    ncol=1,
                    heights = c(2, 1.175),
                    font.label = list(size=9, face="plain"))

ggsave(plot= p1_tot, filename=paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", "SPT_trans.pdf"),
       width=180, height=100, units="mm")
