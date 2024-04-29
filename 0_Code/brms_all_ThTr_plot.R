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
  replace_na(list(Train_TD = 0, Pred_TD = 0,
                  Train_TT=-1, Pred_TT=-1)
  ) %>% select(-Folds,-Case) %>%
  filter(!is.na(SORENSEN))
rm(SORENSEN_model_df)
gc()

TrimThin2_df <- FILTERED_df %>% 
  filter(Subset %in% c("Val"), 
         Pred_TW == Train_TW,
         PREDICTORS %in% c("BV", "MI", "B&I"),
         Pred_TRIM == Train_TRIM,
         Pred_TRIM != "NoTrim",
         Pred_THIN == Train_THIN,
         Pred_TD == Train_TD,
         Pred_TT == Train_TT)  %>%
    mutate(Pred_TD = as.factor(Pred_TD), Pred_TT=as.factor(Pred_TT))

plot_df <- TrimThin2_df %>% group_by(PREDICTORS,Pred_TT, Pred_TD, Pred_TW) %>%
  mutate(SORENSEN=SORENSEN*100) %>%
  summarise(MSø=mean(SORENSEN), SDSø=sd(SORENSEN)) %>% 
  pivot_longer(cols=MSø:SDSø,names_to = "metric", values_to = "SORENSEN") %>%
  mutate(Pred_TT= str_c(Pred_TT, " Y") %>% 
           str_replace_all(., "-1 Y", "spThin"))

TT_levels  <- as.character(c("0 Y","1 Y","5 Y","10 Y","70 Y", "spThin"))

plot_ls <- lapply(1:length(unique(plot_df$metric)), FUN = function(y){
  opts <- c("D", "A")
  titles <- c("MSø", "SDSø")
  x <- unique(plot_df$metric)[y]
  gg<-ggplot(
    data = plot_df[plot_df$metric == x, ],
    aes(x = factor(Pred_TD), 
        y = factor(Pred_TT, levels = TT_levels),
        fill = SORENSEN)) +
    geom_tile() + scale_fill_viridis_c(direction = -1, begin = 0,
                                       breaks = scales::breaks_pretty(n=3),
                                       option = opts[y], name = titles[y]) + 
    theme_light() +
    facet_nested(Pred_TW~metric+PREDICTORS, scales = "free", space="free") +
    theme(
      legend.key.width= unit(6,"mm"),
      legend.key.height = unit(3,"mm"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(color= "black",
                                margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
      legend.title = element_blank(),
      legend.text = element_text(size=8),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.margin = margin(-2, 0, 0, -2, unit = "mm"),
      panel.spacing = unit(1, "mm"),
      legend.position = "bottom",
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
p2 <- ggarrange(plotlist = plot_ls, ncol = 2, widths = c(1,0.85)) #

fname <- "TrimThin_brms.RData"
if(file.exists(file.path(Dir.Outputs, fname))){
  load(file.path(Dir.Outputs, fname))
}else{
  TrimThin2_brms <- brm(SORENSEN ~ PREDICTORS + Pred_TW + 
                          Pred_TD * Pred_TT + 
                          (1|Species)
                        , 
                        data = TrimThin2_df,
                        family = "zero_one_inflated_beta",
                        chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(TrimThin2_brms, file = file.path(Dir.Outputs, fname))
}

pb2<- PlotBRMS2(TrimThin2_brms, "Thin") +
  theme(plot.margin = margin(0,1,0,0, "mm"))+
  scale_x_continuous(breaks = scales::breaks_pretty(n=3))

p2_tot <- ggarrange(plotlist = list(p2, pb2), 
                    labels =c("(a)", "(b)"), 
                    label.x = -0.02,common.legend = FALSE,
                    ncol=2,
                    font.label = list(size=9, face="plain"))


ggsave(p2_tot, filename = file.path(Dir.Outputs, "Preds_all.pdf"),
       width=180, height=80, units="mm")

# ggsave(p2_tot, filename = file.path(Dir.Outputs, "Preds_all_diss.pdf"),
#        width=160, height=80, units="mm")
