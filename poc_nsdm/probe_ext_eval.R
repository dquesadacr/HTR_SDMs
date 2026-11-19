# probe_ext_eval.R — paper Ext (2001-2015) ensemble metrics per predictor
# approach x training config, averaged over years at max_sens_spec. Read-only.
base <- "Full/C1_F10/2_Outputs/0_Model_performance/Ensemble/meanw"
sp <- c("Arnica_montana", "Cicerbita_alpina", "Meum_athamanticum")
suppressMessages(library(dplyr))
for (s in sp) {
  f <- file.path(base, s, "eval_regular_biovars.rds")
  e <- as.data.frame(readRDS(f))
  ex <- e %>% filter(Subset == "Ext", threshold == "max_sens_spec")
  agg <- ex %>% group_by(Predictors, ObsSubset) %>%
    summarise(nyr = n_distinct(Year), AUC = mean(AUC), TSS = mean(TSS),
              SOR = mean(SORENSEN), BOYCE = mean(BOYCE), IMAE = mean(IMAE),
              .groups = "drop") %>% arrange(AUC)
  cat("=====", s, "(", nrow(ex), "Ext rows )\n")
  print(as.data.frame(agg), digits = 3)
}
