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
# setwd("/home/dqc/Documents/PhD/scripts/R/sdm/Collab/sdm_sachsen/")

source("mods.R")
source("aux_proj.R")

# ARGUMENTS TO PARSE =======

Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID"))

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

# READING DATA =================================================================

eval_all <- sapply(list.files(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", Spec4folders),
                              pattern = "*.rds$",  full.names = TRUE) %>% str_subset("eval"),
                   simplify = FALSE, readRDS) %>% bind_rows()

dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")

eval_all_long <- eval_all %>%
  mutate(Train = Train  %>% str_replace_all(dictio_replace),
         Train = Train %>% str_remove_all("Obs"),
         ObsSubset = ObsSubset %>% str_remove_all("Obs"),
         ObsSubset = ObsSubset %>% str_replace_all(dictio_replace),
         Predictors = Predictors %>% str_replace_all(dictio_replace)) %>%
  filter(Train == paste0(str_split_i(Predictors, "_", 1), "_" , ObsSubset)) %>%
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
    EvalThinTemp = fct_relevel(EvalThinTemp %>% str_replace_all("100Y", "70Y"), c("N/A", "0Y", "1Y", "5Y", "10Y", "70Y")), #"1Y",
    TrainThinTemp = fct_relevel(TrainThinTemp %>% str_replace_all("100Y", "70Y"), c("N/A", "0Y", "1Y", "5Y", "10Y", "70Y")), #"1Y",
    EvalThinDist = fct_relevel(EvalThinDist %>% str_replace_all("km", " km"), c("N/A", "1 km", "2 km")),
    TrainThinDist = fct_relevel(TrainThinDist %>% str_replace_all("km", " km"), c("N/A", "1 km", "2 km")),
    Subset = fct_relevel(Subset, c("Train", "Val", "Ext")))
#     Subset = fct_relevel(Subset, c("Train", "Val", "Ext", "Full")))

H2000_subset <- eval_all_long %>%
  filter(metric == "SORENSEN", threshold == "max_sens_spec",
         Posteriori == "Reg", EvalObsApproach == "H2k")


H2000_plot <- ggplot(H2000_subset,
                     aes(x=TrainPredictors, y=value,
                         color=TrainPredictors, shape=Subset)) +
  geom_point(size = 1) +
  scale_shape(name = "") +
  scale_y_continuous(labels=scales::label_number(0.01))+
  facet_nested(EvalTrim~"H2k"+EvalThinDist+EvalThinTemp,
    scales = "free_x",
    nest_line = element_line(color="black", linewidth = .2)) +
  theme_light(base_size = 10) +
  scale_color_jco(name="") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0,0, unit = "mm")),
    legend.position = "bottom",
    legend.text = element_text(margin=margin(0, 0, 0, -1.5, unit = "mm")),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-3.75, 0, -1.5, -4, unit = "mm"),
    legend.box = "vertical",
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(-2, "mm"),
    panel.spacing.x = unit(0.5, "mm"),
    panel.spacing.y = unit(1.35, "mm"),
    plot.margin = margin(0, 0.5, 0, 0, unit = "mm"))

rest_subset <- eval_all_long %>%
  filter(metric=="SORENSEN", threshold =="max_sens_spec", Posteriori=="Reg",
         EvalTrim!="NoTrim",
         TrainPredictors!="WC",TrainPredictors!="WC-R",
         ExPredsApproach!="H2k", ExPredsApproach==TrainApproach
         )

rest_plot <- ggplot(rest_subset ,
                    aes(color=Subset,
                        y=value,
                        x=EvalThinTemp)) +
  geom_boxplot(width = 0.5, outlier.size = 0.4, size=0.4) +
  scale_y_continuous(breaks = c(0.25,0.75), minor_breaks = c(0,0.5,1)) +
  scale_color_jco(name="") +
  facet_nested(
    TrainPredictors~EvalObsApproach+EvalThinDist,
    scales = "free_x", space ="free_x",
    nest_line = element_line(color="black", linewidth = .2)) +
  geom_text(data = rest_subset %>%
              group_by(Subset,EvalObsApproach,TrainPredictors,EvalThin,
                       EvalThinDist,EvalThinTemp) %>%
              summarise(Mean = round(mean(value*100),0)),
            aes(x=EvalThinTemp, y=0.2, label=Mean),
            color = "black",
            position = position_dodge2(width = 1),
            size=2.25) +
  theme_light(base_size = 10) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black",
                              margin = margin(0.5,0,0,0, unit = "mm")),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-3.5, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(0.5, "mm"),
    plot.margin = margin(0, 0, 0, 0, unit = "mm"))

joint_plot <- ggarrange(plotlist = list(H2000_plot, rest_plot),
                        labels =c("(a)", "(b)"), label.x = -0.015,
                        widths = c(1,2.5),
                        font.label = list(size=8, face="plain"))

ggsave(plot = joint_plot,
       filename = paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                         Spec4folders, "/eval.pdf"),
       width=180, height=60, units="mm")
