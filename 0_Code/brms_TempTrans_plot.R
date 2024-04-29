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

source("aux_plots.R")


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
                  , Train_TT=-1, Pred_TT=-1)
  ) %>% select(-Folds,-Case) %>%
  filter(!is.na(SORENSEN))
rm(SORENSEN_model_df)
gc()

TempTrans_df <- FILTERED_df %>% 
  filter(PREDICTORS %in% c("BV", "MI", "B&I"), 
         Subset %in% c("Full", "Val"))

plot_df <- TempTrans_df %>% group_by(PREDICTORS,Pred_TW, Train_TW) %>%
  mutate(SORENSEN=SORENSEN*100) %>%
  summarise(MSø=mean(SORENSEN), SDSø=sd(SORENSEN)) %>% 
  pivot_longer(cols=MSø:SDSø,names_to = "metric", values_to = "SORENSEN") 

plot_ls <- lapply(1:length(unique(plot_df$metric)), FUN = function(y){
  opts <- c("D", "A")
  titles <- c("MSø", "SDSø")
  x <- unique(plot_df$metric)[y]
  gg<-ggplot(
    data = plot_df[plot_df$metric == x, ],
    aes(x = factor(Train_TW, levels = c("H2k", "Dec","R10")),
        y = factor(Pred_TW, levels = c("H2k", "Dec","R10")),
        fill = SORENSEN)) +
    geom_tile() + scale_fill_viridis_c(direction = -1, begin = 0,
                                       breaks = scales::breaks_pretty(n=3),
                                       option = opts[y], name = titles[y]) +
    theme_light() +
    labs(x = "Train", y = "Predict") +
    facet_nested(factor(PREDICTORS, levels = c("WC", "WC-R", "BV", "MI", "B&I"))~metric,
                 scales = "free", space="free") +
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
      plot.margin = margin(0, 0, 0, 0, unit = "mm"))
  if(x=="MSø"){gg <- gg+theme(strip.text.y = element_blank(),
                              plot.margin = margin(0, 0.5, 0, 0, unit = "mm"))
  } else {gg <- gg+theme(axis.text.y = element_blank(),
                         axis.ticks.length.y.left = unit(0, "mm"),
                         axis.title.y = element_blank(),
                         plot.margin = margin(0, 0, 0, 0.5, unit = "mm"))}
  return(gg)
})
p2 <- ggarrange(plotlist = plot_ls, ncol = 2, widths = c(1,0.85)) #

TempTrans_df_fix <- TempTrans_df %>% 
  mutate(Pred_TW = factor(str_c("P-", Pred_TW), 
                               levels = paste0("P-", c("H2k","Dec", "R10"))),
         Train_TW = factor(str_c("T-",Train_TW),
                                levels = paste0("T-", c("H2k","Dec", "R10"))))

fname<-"TempTrans_brms.RData"
if(file.exists(file = file.path(Dir.Outputs, fname))){
  load(file = file.path(Dir.Outputs, fname))
}else{
  TempTrans_brms <- brm(SORENSEN ~ PREDICTORS +
                          Pred_TW * Train_TW +
                          (1|Species)
                        , 
                        data = TempTrans_df_fix,
                        family = "zero_one_inflated_beta",
                        chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(TempTrans_brms, file = file.path(Dir.Outputs, fname))
}

pb2 <- PlotBRMS2(TempTrans_brms) +
  labs(x="Coefficient Estimate") +
  scale_x_continuous(breaks = scales::breaks_pretty(n=3)) +
  theme(plot.margin = margin(0.25,1,0.25,1, "mm")) +
  facet_manual(BLOCK~., design="A\nB\nC",
               scales="free", heights = c(2,2,4),
               strip = strip_nested())

p2_tot <- ggarrange(plotlist = list(p2, pb2), 
                    labels =c("(a)", "(b)"), 
                    label.x = -0.02,common.legend = FALSE, 
                    ncol=2,
                    font.label = list(size=9, face="plain"))

ggsave(p2_tot, filename = file.path(Dir.Outputs, "Temp_trans.pdf"),
       width=180, height=80, units="mm")
# ggsave(p2_tot, filename = file.path(Dir.Outputs, "Temp_trans_diss.pdf"),
#        width=160, height=80, units="mm")
