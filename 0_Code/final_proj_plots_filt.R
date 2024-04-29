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
uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

source("mods.R")
source("aux_proj.R")

# case <- c(cmdArg("c"), cmdArg("case"))
#
# setwd(paste0("/home/dqc/Documents/PhD/papers/GCB/iters8/C",case, "_F10/"))
# setwd(paste0("../C",case, "_F10/"))

dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")
dictio_replace_reverse <- c("N/A"=NA, "70Y"="100Y", "H2k"="H2000", " km"="km",
                            "R10"="run10", "Dec"="dec", "WC$"="worldclim_orig",
                            "WC-R"="worldclim_reproj", "MI"="ind", "B&I"="both", 
                            "BV"="biovars", "^Thin"="ObsThin", "^Trim"="ObsTrim")
dictio_replace_2 <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10", "dec"="Dec")


mss_test <- list.files(pattern = paste0("elev.rds"), recursive = TRUE) %>% 
  str_subset("mascula", negate = TRUE) %>%
  sapply(., simplify = FALSE, FUN=function(x){
    niche_metrics <- readRDS(x)
    niche_metrics$Species <- str_split_i(x, "/", 5) %>% str_replace_all("_", "\n")
    return(niche_metrics)
  }) %>% bind_rows()

mss_filt <- mss_test %>% ungroup %>%
  filter(TrainPredictors %in% c("WC", "BV", "B&I"),
         TrainApproachPeriod != "TrainDec" & PredsPeriod != "PredDec") %>%
  separate_wider_delim(Model, delim = "_", too_few = "align_start", cols_remove = FALSE,
                       names= c("TrainPreds", "TrainApproach", "TrainTrim", 
                                "TrainThin", "TrainThinDist", "TrainThinTemp")) %>%
  replace_na(list(TrainThinTemp="N/A", TrainThinDist="N/A")) %>%
  mutate(
    TrainThinTemp = fct_relevel(TrainThinTemp %>% str_replace_all("100Y", "70Y"), c("N/A", "0Y", "1Y", "5Y", "10Y", "70Y")), #"1Y",
    TrainThinDist = fct_relevel(TrainThinDist %>% str_replace_all("km", " km"), 
                                c("N/A", "1 km", "2 km"))) %>%
  mutate(TrainThinDist = fct_relabel(TrainThinDist, ~ gsub("N/A", "0 km", .x)),
         TrainThinTemp = if_else(TrainThinDist == "0 km", 
                                 if_else(TrainThinTemp=="N/A", 
                                         fct_relabel(TrainThinTemp, ~ gsub("N/A", "", .x)),
                                         TrainThinTemp),
                                 fct_relabel(TrainThinTemp, ~ gsub("N/A", "spThin", .x))))

h2k_df <- mss_filt %>% 
  filter(TrainApproach == "H2k",
         TrainPredictors == "WC",
         TrainThinTemp %in% c("spThin"),
         TrainTrim != "Trim",
         TrainThinDist != "0 km",
  ) %>% mutate(PredsGCM="Obs+SSP585", TrainPredictors="WC-R")

r10_df <- mss_filt %>% 
  filter(TrainApproach == "R10",
         TrainThinTemp == "0Y" & TrainThinDist == "1 km" |
           TrainThinTemp == "1Y" & TrainThinDist == "1 km" |
           TrainThinTemp == "5Y" & TrainThinDist == "1 km" |
           TrainThinTemp == "0Y" & TrainThinDist == "2 km" |
           TrainThinTemp == "1Y" & TrainThinDist == "2 km" |
           TrainThinTemp == "5Y" & TrainThinDist == "2 km",
         PredsPeriod != "PredP20",
         !(PredsPeriod=="PredP20" & PredsGCM =="Train"),
         ) %>%
  mutate(PredsGCM = if_else(PredsGCM=="Train", "ReKIS", "Hist+RCP85"))

mss_fin <- bind_rows(h2k_df, r10_df) %>% 
  mutate(PredsGCM = factor(PredsGCM, c("Obs+SSP585", "ReKIS", "Hist+RCP85")),
         PredsPeriod = str_remove_all(PredsPeriod, "Pred"),
         TrainPredictors=factor(TrainPredictors, c("WC-R", "BV", "B&I")))

mss_summ <- mss_fin %>%
  group_by(TrainPredictors,TrainApproach,PredsPeriod, Species, Year, PredsGCM) %>%
  summarise(mPixels = mean(Pixels), sdPixels = sd(Pixels), 
            mElev = mean(Elev), sdElev = sd(Elev),
            n=n())

sec1 <- help_secondary(mss_summ, primary = mPixels, secondary = mElev,
                       breaks = c(400,700,1000))
eb <- aes(ymin = ifelse(mPixels - 2*sdPixels<0, 0, mPixels - 2*sdPixels), 
          ymax = ifelse(mPixels + 2*sdPixels>1, 1, mPixels + 2*sdPixels))

plot_pixels <- ggplot(mss_summ,
                      aes(x=Year, y=mPixels,
                          color=Species)) +
  geom_line() +
  geom_line(aes(y = sec1$proj(mElev)), linetype="dotted") +
  geom_ribbon(eb, alpha=0.15, size=0.25, show.legend = FALSE) +
  scale_y_continuous(sec.axis = sec1, labels=scales::label_percent(),
                     breaks = c(0,.4,.8)) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  facet_nested(Species~TrainPredictors+PredsGCM,
               scales = "free_x", space = "free_x",
               nest_line = element_line(color="black", linewidth = .25)) +
  theme_light(base_size = 10) +
  scale_color_manual(values = c("#7d5690", "#34acf8", "#6acc68", "#f7b036", "#ef9fef",
                                "#758796", "#1929C8", "#C42503", "#3f822e", "#C1C1C1", "#8F0000"),
                     name="") +
  guides(color=guide_legend(nrow=2,byrow=TRUE, override.aes = list(linewidth=0.75))) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=7, angle = 45, hjust = 1),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0,0, unit = "mm")),
    legend.position = "bottom",
    legend.text = element_text(margin=margin(0, 1, 0, 0, unit = "mm"), face = "italic"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-3, 0, 0.25, -2, unit = "mm"),
    legend.spacing.x = unit(2, "mm"),
    legend.spacing.y = unit(0, "mm"),
    panel.spacing = unit(1, "mm"),
    strip.text.y.right = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "mm"))

ggsave(plot= plot_pixels, 
       filename=paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                       "/sel_fin_2sd_filt.pdf"),
       width=180, height=100, units="mm")

# ggsave(plot= plot_pixels,
#        filename=paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
#                        "/sel_fin_2sd_filt_diss.pdf"),
#        width=160, height=100, units="mm")
