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

source("./aux_plots.R")


# DIRECTORIES ==================================================================
Dir.Base <- getwd()
setwd(Dir.Base)
Dir.Outputs <- file.path(Dir.Base, "3_brms")
mkdirs(Dir.Outputs)

# READING DATA =================================================================
dictio_replace <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10")
dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")

SORENSEN_model_df <- readRDS(file.path(Dir.Base,"SORENSEN_model_df.rds"))

FILTERED_df <- SORENSEN_model_df %>% filter(
  # Case == 4 & Folds == 10 & # I have found best performance at this combination
    threshold=="max_sens_spec",
    Species != "Orchis mascula") %>%
  replace_na(list(Train_TD = 0, Pred_TD = 0
                  , Train_TT=-1, Pred_TT=-1)) %>%
  select(-Folds,-Case) %>%
  filter(!is.na(SORENSEN))

rm(SORENSEN_model_df)
gc()

TrimThin_df <- FILTERED_df %>% 
  filter(Subset %in% c("Val", "Ext"),
         Pred_TW == "H2k",
         Pred_TW == Train_TW,
         Pred_TRIM == Train_TRIM,
         Pred_THIN == Train_THIN,
         Pred_TD == Train_TD,
         Pred_TT == Train_TT,
  ) %>%
  mutate(Train_TD = fct_relevel(as.character(Train_TD), c("0", "1", "2")))

plot_df1 <- TrimThin_df %>% group_by(PREDICTORS, Subset) %>%
  mutate(SORENSEN=SORENSEN*100) %>%
  summarise(MSø=mean(SORENSEN), SDSø=sd(SORENSEN)) %>% 
  pivot_longer(cols=MSø:SDSø,names_to = "metric", values_to = "SORENSEN") %>%
  mutate(PREDICTORS=fct_relevel(PREDICTORS, rev(levels(PREDICTORS))))

plot_df2 <- TrimThin_df %>% 
  filter(Subset %in% c("Val")) %>%
  mutate(SORENSEN=SORENSEN*100) %>%
  mutate(Pred_TD = fct_relevel(str_c(Pred_TD, " km"), c("1 km", "2 km"))) %>% 
  group_by(Pred_TRIM, Pred_THIN, Pred_TD) %>%
  summarise(MSø=mean(SORENSEN), SDSø=sd(SORENSEN)) %>% 
  pivot_longer(cols=MSø:SDSø,names_to = "metric", values_to = "SORENSEN")

lims_both<- bind_rows(plot_df1 %>% ungroup() %>% select(metric, SORENSEN),
                      plot_df2 %>% ungroup() %>% select(metric, SORENSEN)) %>%
  group_by(metric) %>%
  summarise(min=min(SORENSEN), max=max(SORENSEN))

plot_ls1 <- lapply(1:length(unique(plot_df1$metric)), FUN = function(y){
  opts <- c("D", "A")
  titles <- c("MSø", "SDSø")
  x <- unique(plot_df1$metric)[y]
  gg<-ggplot(
    data = plot_df1[plot_df1$metric == x, ],
    aes(x = Subset, 
        y = PREDICTORS, 
        fill = SORENSEN)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1, begin = 0,
                         limits = unlist(lims_both[y,-1]), oob = scales::squish,
                         option = opts[y], name = "", breaks = scales::breaks_pretty(n=3)) +
    theme_light() +
    facet_nested(.~metric, scales = "free", space="free") +
    guides(fill="none")+
    theme(
      legend.key.width= unit(6,"mm"),
      legend.key.height = unit(3,"mm"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(color= "black",
                                margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.margin = margin(-2, 0, 0, -2, unit = "mm"),
      panel.spacing = unit(1, "mm"),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "mm"))
  if(x=="MSø"){gg <- gg+theme(strip.text.y = element_blank(),
                              plot.margin = margin(0, 0.5, 0, 0, unit = "mm"))
  } else {gg <- gg+theme(axis.text.y = element_blank(),
                         axis.ticks.length.y.left = unit(0, "mm"),
                         plot.margin = margin(0, 0, 0, 0.5, unit = "mm"))}
  return(gg)
})

p1 <- ggarrange(plotlist = plot_ls1, ncol = 2, widths = c(1,0.825)) +
  theme(plot.margin = margin(0,0.5,0,0,"mm"))

plot_ls2 <- lapply(1:length(unique(plot_df2$metric)), FUN = function(y){
  opts <- c("D", "A")
  titles <- c("MSø", "SDSø")
  x <- unique(plot_df2$metric)[y]
  gg<-ggplot(
    data = plot_df2[plot_df2$metric == x, ],
    aes(x = Pred_TRIM, 
        y = Pred_THIN, 
        fill = SORENSEN)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1, begin = 0,
                         limits = unlist(lims_both[y,-1]),
                         oob = scales::squish,
                         breaks = scales::breaks_pretty(n=3),
                         option = opts[y], name = "") +
    theme_light() +
    facet_nested(Pred_TD~metric, scales = "free", space="free") +
    theme(
      legend.key.width= unit(6,"mm"),
      legend.key.height = unit(3,"mm"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(color= "black", #size=8,
                                margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.margin = margin(-2, 0, 0, -2, unit = "mm"),
      legend.text = element_text(size = 8), 
      panel.spacing = unit(1, "mm"),
      legend.position = "bottom",
      strip.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "mm"))
  if(x=="MSø"){gg <- gg+theme(strip.text.y = element_blank(),
                              plot.margin = margin(1, 0.5, 0, 0, unit = "mm"))
  } else {gg <- gg+theme(axis.text.y = element_blank(),
                         axis.ticks.length.y.left = unit(0, "mm"),
                         plot.margin = margin(1, 0, 0, 0.5, unit = "mm"))}
  return(gg)
})

p2 <- ggarrange(plotlist = plot_ls2, ncol = 2, widths = c(1,0.825)) +
  theme(plot.margin = margin(0,0.5,0,0,"mm"))

fname <- "H2k_PREDH2k_preds.RData"
if(file.exists(file.path(Dir.Outputs, fname))){
  load(file.path(Dir.Outputs, fname))
}else{
  TrimThin_brms <- brm(SORENSEN ~ Subset + 
                         PREDICTORS + 
                         (1|Species)
                       ,
                       data = TrimThin_df,
                       family = "zero_one_inflated_beta",
                       chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(TrimThin_brms, file = file.path(Dir.Outputs, fname))
}
pb1 <- PlotBRMS2(TrimThin_brms) + 
  theme(axis.title.x = element_blank(),
        plot.margin = margin(0,0.5,0,3,"mm")) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=3))

fname <- "H2k_PREDH2k_tt.RData"
if(file.exists(file.path(Dir.Outputs, fname))){
  load(file.path(Dir.Outputs, fname))
}else{
  TrimThin_brms <- brm(SORENSEN ~ Train_TRIM +
                         Train_THIN * Train_TD +
                         (1|Species)
                       ,
                       data = TrimThin_df %>% filter(Subset=="Val"),
                       family = "zero_one_inflated_beta",
                       chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(TrimThin_brms, file = file.path(Dir.Outputs, fname))
}

pb2 <- PlotBRMS2(TrimThin_brms) +
  theme(plot.margin = margin(0,0.5,0,3,"mm"))+
  scale_x_continuous(breaks = scales::breaks_pretty(n=3))

p1_tot <- ggarrange(plotlist = list(p1, pb1, p2, pb2), 
                    labels =c("(a)", "(b)", "(c)", "(d)"), 
                    label.x = -0.02,common.legend = FALSE,
                    label.y = 1.01,
                    nrow = 2,
                    ncol=2,
                    heights = c(1, 1.1),
                    font.label = list(size=9, face="plain"))

ggsave(p1_tot, filename = file.path(Dir.Outputs, "H2k_all.pdf"),
       width=180, height=80, units="mm")

# ggsave(p1_tot, filename = file.path(Dir.Outputs, "H2k_all_diss.pdf"),
#        width=160, height=80, units="mm")
