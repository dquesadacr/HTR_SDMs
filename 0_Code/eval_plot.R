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

# approach <- "run10" # options run5, run10, run15, H2000 and worldclim (WC)
# preds <- "wc" "bv" # "biovars" # read as BV (bv, biovars) also, and indexes or ind
# Spec_ID <-   "alpina" # 
# Spec_ID <- "ianthus"
# Spec_ID <- "centaurea" #
# Spec_ID <- "ranun" #

# approach <- c(cmdArg("approach"), cmdArg("a")) # "run10" # options run5, run10, run15, H2000 and worldclim (WC)
# preds <- c(cmdArg("p"), cmdArg("preds"), cmdArg("predictors")) # "biovars" # read as BV (bv, biovars) also, and indexes or ind
Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID")) #"centaurea" #

uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

# DIRECTORIES ==================================================================
Dir.Base <- getwd()
setwd(Dir.Base)
#setwd("/home/dqc/Documents/PhD/scripts/R/sdm/Collab/sdm_sachsen/")

# READING DATA =================================================================

eval_all <- sapply(list.files(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", Spec4folders),
                              pattern = "*.rds$",  full.names = TRUE) %>% str_subset("eval"),
                   simplify = FALSE, readRDS) %>% bind_rows()

dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")

eval_all_long <- eval_all %>%
  mutate(Train = Train  %>% str_replace_all(dictio_replace), #) %>% #,
         Train = Train %>% str_remove_all("Obs"),
         ObsSubset = ObsSubset %>% str_remove_all("Obs"),
         ObsSubset = ObsSubset %>% str_replace_all(dictio_replace),
         Predictors = Predictors %>% str_replace_all(dictio_replace))  %>%
  filter(Train == paste0(str_split_i(Predictors, "_", 1), "_" ,ObsSubset)) %>%
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
#   filter(Train == paste0(str_split_i(Predictors, "_", 1), "_" ,ObsSubset)) %>%
  filter(metric=="SORENSEN", threshold =="max_sens_spec", Posteriori == "Reg", 
         # EvalThinTemp != "70Y", TrainThinTemp != "70Y",
         EvalObsApproach=="H2k",)
#          ExPredictors == TrainPredictors, TrainApproach=="H2k")


H2000_plot <- ggplot(H2000_subset, #, Train_preds!="worldclim"
                     aes(x=TrainPredictors, y=value, 
                         color=TrainPredictors, shape=Subset)) +
  geom_point(size = 1) + #varwidth = TRUE width = 0.5, size=0.35
  scale_shape(name="") +
  scale_y_continuous(labels=scales::label_number(0.01))+
  facet_nested(EvalTrim~"H2k"+EvalThinDist+EvalThinTemp, #EvalThin
#   facet_nested(TrainTrim~"H2k"+TrainThinDist+TrainThinTemp, #EvalThin
    scales = "free_x",
    nest_line = element_line(color="black", linewidth = .2)) + 
  theme_light(base_size = 10) +
  scale_color_jco(name="") + #, end=0.9
  # guides(color = guide_legend(nrow=2, byrow = TRUE),
  #        shape = guide_legend(nrow=2, byrow = TRUE)) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0,0, unit = "mm")),
    legend.position = "bottom",
    # legend.key.width = unit(7.5, "mm"),
    # legend.key.height = unit(3, "mm"),
    legend.text = element_text(margin=margin(0, 0, 0, -1.5, unit = "mm")),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-3.75, 0, -1.5, -4, unit = "mm"),
    # legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box = "vertical",
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(-2, "mm"),
    # panel.spacing = unit(0.5, "mm"),
    panel.spacing.x = unit(0.5, "mm"),
    panel.spacing.y = unit(1.35, "mm"),
    # plot.title = element_text(hjust = 0.5, margin = margin(0.5,0,0,0, unit = "mm"), size = 10),
    plot.margin = margin(0, 0.5, 0, 0, unit = "mm"))

rest_subset <- eval_all_long %>%
  filter(metric=="SORENSEN", threshold =="max_sens_spec", Posteriori=="Reg",
         EvalTrim!="NoTrim",
         TrainPredictors!="WC",TrainPredictors!="WC-R",
         ExPredsApproach!="H2k", ExPredsApproach==TrainApproach
         # !(TrainPredictors=="WC" & TrainThinTemp == "N/A") ,!(TrainPredictors=="WC-R" & TrainThinTemp == "N/A"),
         )
# test <- rest_subset %>%
#   group_by(Subset,EvalObsApproach,TrainPredictors,EvalThin,
#            EvalThinDist,EvalThinTemp) %>%
#   summarise(Mean = round(mean(value*100),2))

rest_plot <- ggplot(rest_subset ,
                    aes(color=Subset,
                        y=value,
                        x=EvalThinTemp)) +
  geom_boxplot(width = 0.5, outlier.size = 0.4, size=0.4) + #varwidth = TRUE width = 0.5
  # geom_violin() +
  scale_y_continuous(breaks = c(0.25,0.75), minor_breaks = c(0,0.5,1)) +
  scale_color_jco(name="") + #, end=0.9
  facet_nested(
    TrainPredictors~EvalObsApproach+EvalThinDist, #+EvalThinTemp
    scales = "free_x", space ="free_x",
    nest_line = element_line(color="black", linewidth = .2)) + #, nrow=4 + scales = "free",
  geom_text(data = rest_subset %>%
              group_by(Subset,EvalObsApproach,TrainPredictors,EvalThin, 
                       EvalThinDist,EvalThinTemp) %>%
              summarise(Mean = round(mean(value*100),0)),
            aes(
              # x=interaction(Subset,EvalObsApproach),
              x=EvalThinTemp, 
              y=0.2, label=Mean), 
            # color = rep(c("blue", "black"), 57),
            color = "black",
            # position = position_jitter(width = 0.2, height = 0.1, seed = 3),
            position = position_dodge2(width = 1),
            size=2.25) +
  theme_light(base_size = 10) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(size=7), #angle = 45, hjust = 1, 
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", 
                              margin = margin(0.5,0,0,0, unit = "mm")),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-3.5, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(0.5, "mm"),
    # plot.title = element_text(hjust = 0.5, margin = margin(1.5,0,0,0, unit = "mm"), size = 10),
    plot.margin = margin(0, 0, 0, 0, unit = "mm"))


joint_plot <- ggarrange(plotlist = list(H2000_plot, rest_plot),
                        labels =c("(a)", "(b)"), label.x = -0.015,
                        widths = c(1,2.5), #2.75
                        # align = "h",
                        font.label = list(size=8, face="plain"))

ggsave(plot= joint_plot, 
       filename=paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                       Spec4folders, "/eval.pdf"), 
       width=180, height=60, units="mm")

# best <-  eval_all_long %>% 
#   ungroup() %>%
#   filter(metric=="SORENSEN", threshold =="max_sens_spec", 
#          Posteriori == "Reg", TrainPredictors != "WC", 
#          TrainPredictors != "WC-R") %>%
#   # select(!c(n_presences, n_absences, thr_value)) %>%
#   # group_by(EvalObsApproach,TrainPredictors, EvalTrim,
#   #          EvalThin, EvalThinDist,EvalThinTemp) %>%
#   # filter()
#   pivot_wider(names_from = Subset, values_from = c(n_presences, n_absences, thr_value, value)) %>%
#   rowwise %>%
#   mutate(value_Train = case_when(!is.na(value_Val) | !is.na(value_Ext) ~ NA_real_,
#                                  TRUE ~ value_Train),
#          # value = rowMeans(select(., c(value_Train, value_Val, value_Ext)), na.rm=TRUE), 
#          value = mean(c_across(c(starts_with('value'))), na.rm=TRUE),
#          # value = case_when(is.na(value) ~ value_Train,TRUE ~ value)
#          ) #%>%
# 
# 
# best <-  eval_all_long %>% 
#   ungroup() %>%
#   filter(metric=="SORENSEN", threshold =="max_sens_spec", 
#          Posteriori == "Reg", TrainPredictors != "WC", TrainPredictors != "WC-R") %>%
#   select(!c(n_presences, n_absences, thr_value)) %>%
#   # group_by(EvalObsApproach,TrainPredictors, EvalTrim,
#   #          EvalThin, EvalThinDist,EvalThinTemp) %>%
#   # filter()
#   # pivot_wider(names_from = Subset, values_from = c(n_presences, n_absences, thr_value, value)) %>%
#   pivot_wider(names_from = Subset, values_from = c(value)) %>%
#   # rowwise %>%
#   # mutate(value_Train = case_when(!is.na(value_Val) | !is.na(value_Ext) ~ NA_real_,
#   #                                TRUE ~ value_Train),
#   #        # value = rowMeans(select(., c(value_Train, value_Val, value_Ext)), na.rm=TRUE), 
#   #        # value = mean(c_across(c(starts_with('value'))), na.rm=TRUE),
#   #        # value = case_when(is.na(value) ~ value_Train,TRUE ~ value)
#   # ) #%>%
#   rowwise %>%
#   mutate(Train = case_when(!is.na(Val) | !is.na(Ext) ~ NA_real_,
#                                  TRUE ~ Train),
#          value = mean(c(Train, Val, Ext), na.rm=TRUE)) %>%
#   group_by(EvalObsApproach,TrainPredictors, EvalTrim,
#            EvalThin, EvalThinDist,EvalThinTemp) %>%
#   summarise(Mean = round(mean(value*100, na.rm = TRUE),2))
#          
# 
# test1 <- best %>% filter(is.na(value_Train))
# test2 <- best %>% filter(is.na(value_Train))
# 
#   # filter(is.na(value_Train) & !is.na(value_Val))# %>%
#   # ungroup() %>%
#   # drop_na(value_Train)# %>%
#   group_by(EvalObsApproach,TrainPredictors, EvalTrim,
#            EvalThin, EvalThinDist,EvalThinTemp) %>% nest
#   # filter(!is.na(value_Val) & !is.na(value_Ext))
#   summarise(Mean = round(mean(value*100),2)) %>%
# # , Mean_thr = round(mean(thr_value*100),2)) #%>%
#   group_by(TrainPredictors, EvalObsApproach) %>%
#   arrange(desc(Mean)) %>%
#   slice_head(n=2)
# 
# best <- bind_rows(best %>% filter(EvalObsApproach == "H2k") %>% group_by(EvalObsApproach, EvalTrim, EvalThin) %>% slice_head(n=1),
#                   best %>% filter(EvalObsApproach != "H2k") %>% group_by(EvalObsApproach, EvalThinDist) %>% slice_head(n=2))
# 
#   group_by(EvalObsApproach, EvalTrim, EvalThin) %>%
#   # group_by(EvalObsApproach, EvalTrim, EvalThinDist, EvalThinTemp) %>%
#   slice_head(n=2)
#   # summarise(Mean2=mean(Mean))

# best <-  eval_all_long %>%
#   filter(metric=="SORENSEN", threshold =="max_sens_spec",
#          Posteriori == "Reg", TrainPredictors != "WC", TrainPredictors != "WC-R") %>%
#   group_by(Subset,EvalObsApproach,TrainPredictors,EvalThin,
#            EvalTrim,EvalThinDist,EvalThinTemp) %>%
#   summarise(Mean = round(mean(value*100),2)) %>%
#   pivot_wider(names_from = Subset, values_from = c(Mean)) %>%
#   ungroup() %>%
#   rowwise() %>%
#   mutate(Train = case_when(!is.na(Val) | !is.na(Ext) ~ NA_real_,
#                            TRUE ~ Train),
#          Mean = mean(c(Train, Val, Ext), na.rm=TRUE)) %>%
#   group_by(TrainPredictors, EvalObsApproach, EvalThin) %>%
#   arrange(desc(Mean)) # %>%
#   #slice_head(n=2)
#
# best_fin <- bind_rows(best %>% filter(EvalObsApproach == "H2k") %>% group_by(EvalTrim, EvalThin, EvalThinTemp) %>% slice_head(n=1),
#                   best %>% filter(EvalObsApproach == "Dec") %>% group_by(EvalObsApproach, EvalThinDist) %>% slice_head(n=2),
#                   best %>% filter(EvalObsApproach == "R10", EvalThinDist=="N/A") %>% group_by(EvalObsApproach) %>% slice_head(n=1),
#                   best %>% filter(EvalObsApproach == "R10") %>% group_by(EvalObsApproach, TrainPredictors) %>% slice_head(n=2))


# test10 <- best %>%
#   group_by(EvalObsApproach,TrainPredictors, EvalTrim,
#            EvalThin, EvalThinDist,EvalThinTemp) %>%
#   summarise(Mean = round(mean(value*100),2)) %>%
#   # , Mean_thr = round(mean(thr_value*100),2)) #%>%
#   # group_by(EvalObsApproach, EvalTrim, EvalThin) %>%
#   group_by(TrainPredictors, EvalObsApproach) %>%
#   arrange(desc(Mean)) %>% 
#   slice_head(n=2)
# # summarise(Mean2=mean(Mean))
#   
  
# saveRDS(best_fin, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
#                      Spec4folders, "/best_models.rds"))
