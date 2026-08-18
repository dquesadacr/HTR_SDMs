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
  threshold=="max_sens_spec",
  Species != "Orchis mascula") %>%
  replace_na(list(Train_TD = 0, Pred_TD = 0,
                  Train_TT=-1, Pred_TT=-1)
  ) %>% select(-Case) %>%
  filter(!is.na(SORENSEN))

# DEBUG: Check data structure and counts
cat("\n=== DEBUG INFO ===\n")
cat("Total rows in FILTERED_df:", nrow(FILTERED_df), "\n")
cat("Unique Species:", length(unique(FILTERED_df$Species)), "\n")
cat("Unique Folds:", length(unique(FILTERED_df$Folds)), "\n")
cat("SORENSEN range:", min(FILTERED_df$SORENSEN, na.rm=TRUE), "to", max(FILTERED_df$SORENSEN, na.rm=TRUE), "\n\n")

rm(SORENSEN_model_df)
gc()


# =============================================================================
# PANEL A: H2k Analysis - PREDICTORS by Subset (from brms_H2k_plot.R p1 part 1)
# =============================================================================

H2k_df <- FILTERED_df %>% 
  filter(Subset %in% c("Val", "Ext"),
         Pred_TW == "H2k",
         Pred_TW == Train_TW,
         Pred_TRIM == Train_TRIM,
         Pred_THIN == Train_THIN,
         Pred_TD == Train_TD,
         Pred_TT == Train_TT) %>%
  mutate(Train_TD = fct_relevel(as.character(Train_TD), c("0", "1", "2")))

cat("H2k_df rows:", nrow(H2k_df), "\n")
cat("H2k_df by PREDICTORS x Subset:\n")
print(table(H2k_df$PREDICTORS, H2k_df$Subset))

# Boxplot for H2k: PREDICTORS by Subset (using raw SORENSEN values)
# Note: Data has multiple species per combination, so boxplots will show variation across species
H2k_plot_df1 <- H2k_df %>% 
  mutate(SORENSEN_pct = SORENSEN*100)

lims_H2k <- H2k_plot_df1 %>% 
  summarise(min=min(SORENSEN_pct), max=max(SORENSEN_pct)) %>%
  pull()

q_range <- c(0.05, 0.95)
q_range <- c(0.1, 0.9) # Adjusted for better visualization
qs <- quantile(H2k_plot_df1$SORENSEN_pct, q_range, na.rm=TRUE)

cat("\nH2k_plot_df1 rows:", nrow(H2k_plot_df1), "\n")
cat("H2k_plot_df1 by PREDICTORS x Subset:\n")
print(table(H2k_plot_df1$PREDICTORS, H2k_plot_df1$Subset))
cat("\nSORENSEN_pct range:", min(H2k_plot_df1$SORENSEN_pct), "to", max(H2k_plot_df1$SORENSEN_pct), "\n")

p_H2k_preds <- ggplot(data = H2k_plot_df1, 
                      aes(x = PREDICTORS, y = SORENSEN_pct, color = PREDICTORS)) +
  geom_boxplot(outlier.alpha = 0.2) +
  # geom_violin(alpha = 0, draw_quantiles = c(0.25, 0.5, 0.75)) +
  # scale_color_viridis_d(option = "H", name = "", begin = 0.075, end = 0.925) + 
  scale_color_manual(values = manual_colors, name= "") +
  # scale_color_viridis_d(option = "D", name = "", begin = 0, end = 0.75) +  
  theme_light() +
  facet_nested(""~Subset, scales = "free_x", space="free_x") +
  # labs(x = "", y = "Sørensen (%)") +
  labs(y="") +
  theme(
    legend.key.width= unit(6,"mm"),
    legend.key.height = unit(4,"mm"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
    # legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-4, 0, 0, 0, unit = "mm"),
    # legend.box.margin = margin(0, 0, 0, -2, unit = "mm"),
    panel.spacing = unit(2, "mm"),
    legend.position = "bottom", # "top"
    legend.text = element_text(size=8, face = "plain"),
    # legend.title = element_text(size=8, face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.title.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "mm")) +
  scale_y_continuous(breaks = scales::breaks_pretty(n=3), limits = qs) # +
  # scale_y_continuous(limits = lims_H2k)

# Boxplot for H2k: TRIM x THIN by TD (p2 from original script)
H2k_plot_df2 <- H2k_df %>% 
  filter(Subset %in% c("Val")) %>%
  mutate(SORENSEN_pct = SORENSEN*100) %>%
  mutate(Pred_TD = fct_relevel(str_c(Pred_TD, " km"), c("1 km", "2 km"))) %>%
  mutate(Pred_THIN = as.factor(Pred_THIN))

qs2 <- quantile(H2k_plot_df2$SORENSEN_pct, q_range, na.rm=TRUE)

cat("\nH2k_plot_df2 rows:", nrow(H2k_plot_df2), "\n")
cat("H2k_plot_df2 by Pred_TRIM x Pred_THIN:\n")
print(table(H2k_plot_df2$Pred_TRIM, H2k_plot_df2$Pred_THIN))

p_H2k_tt <- ggplot(data = H2k_plot_df2, 
                   aes(x = Pred_THIN, y = SORENSEN_pct, color = Pred_THIN)) +
  geom_boxplot(outlier.alpha = 0.2) +
  # scale_color_viridis_d(option = "H", name = "", begin = 0.075, end = 0.925) + 
  scale_color_manual(values = manual_colors[c(4,5)], name= "") +
  theme_light() +
  facet_nested(Pred_TD~Pred_TRIM, 
               nest_line = element_line(color="black", linewidth = 0.25)) +
              #  scales = "free_x", space="free_x") +
  # labs(y = "Sørensen (%)") + # x = ""
  labs(y="") +
  theme(
    legend.key.width= unit(6,"mm"),
    legend.key.height = unit(4,"mm"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
    # legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-4, 0, 0, 0, unit = "mm"),
    legend.text = element_text(size=8, face = "plain"),
    legend.title = element_text(size=8, face = "plain"),
    panel.spacing = unit(2, "mm"),
    legend.position = "bottom",
    # strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),    
    plot.margin = margin(0, 0, 0, 0, unit = "mm")) +
  scale_y_continuous(breaks = scales::breaks_pretty(n=3), limits = qs2) #+
  # scale_y_continuous(limits = lims_H2k)

# H2k BRMS models
fname_H2k1 <- "H2k_PREDH2k_preds.RData"
if(file.exists(file.path(Dir.Outputs, fname_H2k1))){
  load(file.path(Dir.Outputs, fname_H2k1))
  # The loaded object from original script is named 'TrimThin_brms', assign to our variable
  H2k_brms_preds <- TrimThin_brms
}else{
  H2k_brms_preds <- brm(SORENSEN ~ Subset + 
                          PREDICTORS + 
                          (1|Species)
                        ,
                        data = H2k_df,
                        family = "zero_one_inflated_beta",
                        chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(H2k_brms_preds, file = file.path(Dir.Outputs, fname_H2k1))
}

fname_H2k2 <- "H2k_PREDH2k_tt.RData"
if(file.exists(file.path(Dir.Outputs, fname_H2k2))){
  load(file.path(Dir.Outputs, fname_H2k2))
  # The loaded object from original script is named 'TrimThin_brms', assign to our variable
  H2k_brms_tt <- TrimThin_brms
}else{
  H2k_brms_tt <- brm(SORENSEN ~ Train_TRIM +
                       Train_THIN * Train_TD +
                       (1|Species)
                     ,
                     data = H2k_df %>% filter(Subset=="Val"),
                     family = "zero_one_inflated_beta",
                     chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(H2k_brms_tt, file = file.path(Dir.Outputs, fname_H2k2))
}


# =============================================================================
# TEMP TRANS ANALYSIS (from brms_TempTrans_plot.R) - will be panels c/d
# =============================================================================

TempTrans_df <- FILTERED_df %>% 
  filter(PREDICTORS %in% c("BV", "MI", "B&I"), 
         Subset %in% c("Full", "Val"))

cat("\nTempTrans_df rows:", nrow(TempTrans_df), "\n")
cat("TempTrans_df by PREDICTORS x Train_TW x Pred_TW:\n")
print(table(TempTrans_df$PREDICTORS, TempTrans_df$Train_TW, TempTrans_df$Pred_TW))

TempTrans_plot_df <- TempTrans_df %>% 
  mutate(SORENSEN_pct = SORENSEN*100)

qs3 <- quantile(TempTrans_plot_df$SORENSEN_pct, q_range, na.rm=TRUE)

p_TempTrans <- ggplot(
  data = TempTrans_plot_df,
  aes(x = Train_TW, y = SORENSEN_pct, color = Train_TW)) +
  geom_boxplot(outlier.alpha = 0.2) +
  # scale_color_viridis_d(option = "H", name = "", begin = 0.075, end = 0.925) + 
  scale_color_manual(values = manual_colors, name= "") +
  theme_light() +
  # add legend title

  labs(x = "", y = "") + # , y = "Sørensen (%)"
  guides(color = guide_legend(title = "Train", )) +
  facet_nested(factor(PREDICTORS, levels = c("WC", "WC-R", "BV", "MI", "B&I"))~"Predict"+Pred_TW,
               scales = "free_x", space="free_x",
               nest_line = element_line(color="black", linewidth = 0.25)) +
  theme(
    legend.key.width= unit(6,"mm"),
    legend.key.height = unit(4,"mm"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
    # legend.title = element_blank(),
    legend.text = element_text(size=8, face = "plain"),
    legend.title = element_text(size=8, face = "plain"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-4, 0, 0, 0, unit = "mm"),
    # axis.title.x = element_text(size = 8, color="black",
    #                             margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
    axis.title.y = element_text(size = 8),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),    
    panel.spacing = unit(2, "mm"),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0, unit = "mm")) +
  scale_y_continuous(breaks = scales::breaks_pretty(n=2), limits = qs3)

# TempTrans BRMS model
TempTrans_df_fix <- TempTrans_df %>% 
  mutate(Pred_TW = factor(str_c("P-", Pred_TW), 
                          levels = paste0("P-", c("H2k","Dec", "R10"))),
         Train_TW = factor(str_c("T-",Train_TW),
                           levels = paste0("T-", c("H2k","Dec", "R10"))))

fname_TempTrans <- "TempTrans_brms.RData"
if(file.exists(file = file.path(Dir.Outputs, fname_TempTrans))){
  load(file = file.path(Dir.Outputs, fname_TempTrans))
  # The loaded object from original script is named 'TempTrans_brms', assign to our variable
  TempTrans_brms <- TempTrans_brms
}else{
  TempTrans_brms <- brm(SORENSEN ~ PREDICTORS +
                          Pred_TW * Train_TW +
                          (1|Species)
                        , 
                        data = TempTrans_df_fix,
                        family = "zero_one_inflated_beta",
                        chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(TempTrans_brms, file = file.path(Dir.Outputs, fname_TempTrans))
}


# =============================================================================
# THIN-THIN ANALYSIS (from brms_all_ThTr_plot.R) - will be panels e/f
# =============================================================================

TT_levels <- as.character(c("spThin", "0 Y","1 Y","5 Y","10 Y","70 Y"))

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

cat("\nTrimThin2_df rows:", nrow(TrimThin2_df), "\n")
cat("TrimThin2_df by PREDICTORS x Pred_TW:\n")
print(table(TrimThin2_df$PREDICTORS, TrimThin2_df$Pred_TW))


TrimThin2_plot_df <- TrimThin2_df %>% 
  mutate(SORENSEN_pct = SORENSEN*100) %>%
  mutate(Pred_TT= str_c(Pred_TT, " Y") %>% 
           str_replace_all(., "-1 Y", "spThin")) %>%
  mutate(Pred_TD = fct_relevel(str_c(Pred_TD, " km"), c("1 km", "2 km")))

qs4 <- quantile(TrimThin2_plot_df$SORENSEN_pct, q_range, na.rm=TRUE)

# Debug: Check for NA values
cat("ThTr SORENSEN summary:", summary(TrimThin2_plot_df$SORENSEN), "\n")
cat("ThTr SORENSEN_pct summary:", summary(TrimThin2_plot_df$SORENSEN_pct), "\n")

# Get limits for ThTr first - check data before calculating limits
lims_ThTr <- c(
  min(TrimThin2_plot_df$SORENSEN_pct, na.rm=TRUE),
  max(TrimThin2_plot_df$SORENSEN_pct, na.rm=TRUE)
)

cat("ThTr limits:", lims_ThTr[1], "to", lims_ThTr[2], "\n")

# If limits are invalid, use default range
if(is.na(lims_ThTr[1]) || is.na(lims_ThTr[2]) || lims_ThTr[1] >= lims_ThTr[2]){
  cat("WARNING: Invalid ThTr limits detected, using default range 0-100\n")
  lims_ThTr <- c(0, 100)
}

p_ThTr <- ggplot(
  data = TrimThin2_plot_df,
  aes(x = factor(Pred_TT, levels = TT_levels), 
      y = SORENSEN_pct, color = Pred_TD)) +
  geom_boxplot(outlier.alpha = 0.2, position = position_dodge(width = 0.95)) +
  # scale_color_viridis_d(option = "H", name = "", begin = 0.075, end = 0.925) + 
  scale_color_manual(values = manual_colors[c(4,5)], name= "") +
  scale_x_discrete(drop=TRUE) +
  # scale_y_continuous(limits = lims_ThTr) +
  theme_light() +
  facet_nested(PREDICTORS~Pred_TW, scales = "free_x", space="free_x",
               nest_line = element_line(color="black", linewidth = 0.25)) +
  labs(x = "", y = "") + # Sørensen (%)
  theme(
    legend.key.width= unit(6,"mm"),
    legend.key.height = unit(4,"mm"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
    legend.title = element_blank(),
    legend.text = element_text(size=8, face = "plain"),
    # legend.title = element_text(size=8, face = "plain"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-4, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(2, "mm"),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 37.5, hjust = 1),    
    plot.margin = margin(0, 0, 0, 0, unit = "mm")) +
  scale_y_continuous(breaks = scales::breaks_pretty(n=3), limits = qs4)

# ThTr BRMS model
fname_ThTr <- "TrimThin_brms.RData"
if(file.exists(file.path(Dir.Outputs, fname_ThTr))){
  load(file.path(Dir.Outputs, fname_ThTr))
  # The loaded object from original script is named 'TrimThin2_brms', assign to our variable
  TrimThin2_brms <- TrimThin2_brms
}else{
  TrimThin2_brms <- brm(SORENSEN ~ PREDICTORS + Pred_TW + 
                          Pred_TD * Pred_TT + 
                          (1|Species)
                        , 
                        data = TrimThin2_df,
                        family = "zero_one_inflated_beta",
                        chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(TrimThin2_brms, file = file.path(Dir.Outputs, fname_ThTr))
}


# R10 figure

r10_df <- FILTERED_df %>% 
  filter(PREDICTORS %in% c("BV", "MI", "B&I"),
         Subset == "Val",
         Pred_TW == Train_TW,
         Pred_TRIM == Train_TRIM,
         Pred_THIN == Train_THIN,
         Train_THIN == "Thin",
         Pred_TD == Train_TD,
         Pred_TT == Train_TT,
         Train_TW == "R10") %>%
  mutate(Pred_TD = as.factor(Pred_TD), Pred_TT=as.factor(Pred_TT))

r10_plot_df <- r10_df %>%
  mutate(SORENSEN_pct = SORENSEN*100) %>%
  mutate(Pred_TT= factor( str_c(Pred_TT, " Y") %>%
         str_replace_all(., "-1 Y", "spThin"), levels =TT_levels)) %>%
  mutate(Pred_TD = fct_relevel(str_c(Pred_TD, " km"), c("1 km", "2 km")))

qs5 <- quantile(r10_plot_df$SORENSEN_pct, q_range, na.rm=TRUE)


p_r10 <- ggplot(
  data = r10_plot_df,
  aes(x = Pred_TT, 
      y = SORENSEN_pct, color = Pred_TT)) +
  geom_boxplot(outlier.alpha = 0.2, position = position_dodge(width = 0.95)) +
  scale_color_manual(values = manual_colors, name= "") + #[c(4,5)]
  scale_x_discrete(drop=TRUE) +
  theme_light() +
  facet_nested(Pred_TD~PREDICTORS, scales = "free_x", space="free_x",
               nest_line = element_line(color="black", linewidth = 0.25)) +
  labs(x = "", y = "") + # Sørensen (%)
  theme(
    legend.key.width= unit(6,"mm"),
    legend.key.height = unit(4,"mm"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", size=8,
                              margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
    legend.title = element_blank(),
    legend.text = element_text(size=8, face = "plain"),
    # legend.title = element_text(size=8, face = "plain"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(-4, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(2, "mm"),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(0, 0, 0, 0, unit = "mm")) +
  scale_y_continuous(breaks = scales::breaks_pretty(n=3), limits = qs5)

# ggsave(p_r10,
#        filename = file.path(Dir.Outputs, "r10_test.pdf"),
#        width=180, height=100, units="mm")

# ThTr BRMS model
fname_r10 <- "R10_brms.RData"
if(file.exists(file.path(Dir.Outputs, fname_r10))){
  load(file.path(Dir.Outputs, fname_r10))
  # The loaded object from original script is named 'R10_brms', assign to our variable
  r10_brms <- r10_brms
}else{
  r10_brms <- brm(SORENSEN ~ PREDICTORS + Pred_TW + 
                    Pred_TD * Pred_TT + 
                    (1|Species)
                  , 
                  data = R10_df,
                  family = "zero_one_inflated_beta",
                  chains = 6, cores = 6, thin = 2, iter = 7e3)
  save(r10_brms, file = file.path(Dir.Outputs, fname_r10))
}


# =============================================================================
# COMBINED FIGURE - 2 COLUMN LAYOUT (pairs side-by-side)
# =============================================================================

# Create BRMS plots
pb_H2k_preds <- PlotBRMS2(H2k_brms_preds, xlabel=FALSE, intercept_text="Val\nWC") +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,0,0,3,"mm")) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=3),
                     labels = prettyZero)

pb_H2k_tt <- PlotBRMS2(H2k_brms_tt, xlabel=FALSE, intercept_text="NoTrim\n1 km\nspThin") +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,0,0,3,"mm"))+
  scale_x_continuous(breaks = scales::breaks_pretty(n=3),
                     labels = prettyZero)

pb_ThTr <- PlotBRMS2(TrimThin2_brms, "Thin", intercept_text="H2k\nBV\nspThin\n1 km") +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
  plot.margin = margin(0,0,0,3, "mm"))+
  scale_x_continuous(breaks = scales::breaks_pretty(n=3),
                     labels = prettyZero) +
  facet_manual(BLOCK~., design="AB\nCD",
               scales="free", heights = c(1,2.85), #4
               widths = c(1,1), #1
               strip = strip_nested())

pb_TempTrans <- PlotBRMS2(TempTrans_brms, intercept_text="PrH2k\nBV\nTrH2k") +
  theme(#axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,0,0,3,"mm")) + #,
  # labs(x="Coefficient Estimate") +
  scale_x_continuous(breaks = scales::breaks_pretty(n=3),
                     labels = prettyZero) 
  #                    +
  # facet_manual(BLOCK~., 
  #             #  design="A\nB\nC",
  #             #  scales="free", heights = c(2,2,7), #4
  #              design="AB\nCC",
  #              scales="free", heights = c(2,3), widths = c(1,1), #1
  #              strip = strip_nested())

pb_r10 <- PlotBRMS2(r10_brms, intercept_text="BV\n1 km\n0 Y") +
  theme( #axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,0,0,3,"mm"))+
  scale_x_continuous(breaks = scales::breaks_pretty(n=3),
                     labels = prettyZero) #+
  # facet_manual(BLOCK~., design="AB\nCC",
  #              scales="free", heights = c(2,3), #4
  #              widths = c(1,1), #1
  #              strip = strip_nested())               

# Combine into pairs: boxplot on left (a,c,e,g), BRMS on right (b,d,f,h)

p_row1 <- ggarrange(
  plotlist = list(p_H2k_preds + 
                  theme(plot.margin = margin(0,2,0,0, "mm"),
                        panel.grid.minor = element_blank()),
                  pb_H2k_preds),
  labels = c("(a)", "(b)"),
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 2,
  widths = c(1, 1), #.1
  font.label = list(size=9, face="plain")
)

p_row2 <- ggarrange(
  plotlist = list(p_H2k_tt + 
                  theme(plot.margin = margin(0,2,0,0, "mm"),
                        panel.grid.minor = element_blank()),
                  pb_H2k_tt),
  labels = c("(c)", "(d)"),
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 2,
  widths = c(1, 1), #.1
  font.label = list(size=9, face="plain")
)

p_row3 <- ggarrange(
  plotlist = list(p_ThTr + 
                  theme(plot.margin = margin(0,2,0,0, "mm"),
                        panel.grid.minor = element_blank()),
                  pb_ThTr),
  labels = c("(e)", "(f)"),
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 2,
  widths = c(1, 1), #.1
  font.label = list(size=9, face="plain")
)

p_row4 <- ggarrange(
  plotlist = list(p_r10 + 
                  theme(plot.margin = margin(0,2,0,0, "mm"),
                        panel.grid.minor = element_blank()),
                  pb_r10 + theme(axis.title.x = element_blank())),
  labels = c("(g)", "(h)"),
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 2,
  widths = c(1, 1), #.1
  font.label = list(size=9, face="plain")
)

p_row5 <- ggarrange(
  plotlist = list(p_TempTrans + 
                  theme(plot.margin = margin(0,2,0,0, "mm"),
                        panel.grid.minor = element_blank()),
                  pb_TempTrans),
  # labels = c("(g)", "(h)"),
  labels = c(" (i)", "(j)"),
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 2,
  widths = c(1, 1), #.1
  font.label = list(size=9, face="plain")
)

p_row5b <- ggarrange(
  plotlist = list(p_TempTrans + 
                  theme(plot.margin = margin(0,2,0,0, "mm"),
                        panel.grid.minor = element_blank()),
                  pb_TempTrans),
  labels = c("(g)", "(h)"),
  # labels = c(" (i)", "(j)"),
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 2,
  widths = c(1, 1), #.1
  font.label = list(size=9, face="plain")
)

# Final combined figure: all 4 rows stacked vertically (2 columns each)
p_final <- ggarrange(
  plotlist = list(p_row1 + theme(plot.margin = margin(1,1,2,0, "mm")), 
                  p_row2 + theme(plot.margin = margin(0,1,2,0, "mm")), 
                  p_row3 + theme(plot.margin = margin(0,1,2,0, "mm")), 
                  p_row4  + theme(plot.margin = margin(0,1,2,0, "mm")),
                  p_row5  + theme(plot.margin = margin(0,1,0,0, "mm"))
                  ),
  # labels = c("A", "B", "C", "D"),
  # heights = c(1.2, 1.9, 3.1, 2.8),
  # nrow = 4,
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 1,
  # heights = c(1, 1.25, 2.25, 2.25),
  nrow = 5,
  # labels = c("A", "B", "C", "D", "E"),
  heights = c(1.2, 1.9, 3.1, 1.9, 2.7),
  font.label = list(size=9, face="plain") #,
  # spacing = 1
)

# Save the combined figure
ggsave(p_final, 
       filename = file.path(Dir.Outputs, "brms_all_combined+r10.pdf"),
       width=180, 
      #  height=180,
       height=200,
       units="mm")


cat("Combined figure saved to:", file.path(Dir.Outputs, "brms_all_combined.pdf"), "\n")

p_final2 <- ggarrange(
  plotlist = list(p_row1 + theme(plot.margin = margin(1,1,2,0, "mm")), 
                  p_row2 + theme(plot.margin = margin(0,1,2,0, "mm")), 
                  p_row3 + theme(plot.margin = margin(0,1,2,0, "mm")), 
                  # p_row4  + theme(plot.margin = margin(0,1,2,0, "mm")),
                  p_row5b  + theme(plot.margin = margin(0,1,0,0, "mm"))
                  ),
  # labels = c("A", "B", "C", "D"),
  heights = c(1.2, 1.9, 3.1, 2.7),
  nrow = 4,
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 1,
  font.label = list(size=9, face="plain") #,
  # spacing = 1
)

ggsave(p_final2, 
       filename = file.path(Dir.Outputs, "brms_all_combined.pdf"),
       width=180, 
      #  height=180,
       height=180,
       units="mm")

r10_only <- ggarrange(
  plotlist = list(p_r10 + 
                  theme(plot.margin = margin(0,2,0,0, "mm"),
                        panel.grid.minor = element_blank()),
                  pb_r10),
  labels = c("(a)", "(b)"),
  # label.x = -0.02,
  hjust = 0,
  common.legend = FALSE,
  ncol = 2,
  widths = c(1, 1), #.1
  font.label = list(size=9, face="plain")
)

ggsave(r10_only, 
       filename = file.path(Dir.Outputs, "R10.pdf"),
       width=180, height=60, units="mm")