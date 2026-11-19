#!/usr/bin/env Rscript
# 07_stage2_figs.R — stage 2 figures in the paper's own visual language:
# theme_light(base_size=10), ggsci jco palette (scale_color_jco), nested facets,
# turbo viridis fills, 180-mm ggsave, micro-margins — copied from eval_plot.R
# and thin_plot.R, not invented.
#
# fig: stage2_eval.pdf  — Ext metrics (max_sens_spec): x = model, color = model,
#      shape = species, nested facet by metric group (binary / ranking).
# fig: stage2_maps.pdf  — Cicerbita suitability surfaces A | B | C3w | C4 +
#      occurrence layer, turbo fill, per the paper's map style.
#
# Usage: Rscript 07_stage2_figs.R <POC_DIR>
suppressMessages({library(tidyverse); library(terra); library(sf)
                  library(tidyterra); library(ggsci); library(ggh4x)
                  library(ggpubr); library(patchwork)})
POC <- commandArgs(TRUE)[1]
OUT <- paste0(POC, "/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dictio <- c("A_paper" = "A (local HTR)", "B_broad" = "B (whole-range WC)",
            "C3w" = "C3w (meanw N-SDM)", "C3s" = "C3s (avg)", "C4" = "C4 (B covariate)")
grp <- c("A (local HTR)" = "binary", "B (whole-range WC)" = "binary",
         "C3w (meanw N-SDM)" = "binary", "C3s (avg)" = "binary",
         "C4 (B covariate)" = "binary")

## ---- fig a: metric panel ----------------------------------------------------
files <- list.files(POC, "stage2_eval_ext[.]csv$", full.names = TRUE, recursive = TRUE)
stopifnot(length(files) >= 1)
ev <- files %>% map_dfr(read_csv, show_col_types = FALSE) %>%
  mutate(species = str_replace_all(species, "_", " "),
         model = factor(model %>% recode(!!!dictio),
                        names(dictio)),
         species = factor(species, c("Arnica montana", "Cicerbita alpina",
                                     "Meum athamanticum"))) %>%
  filter(threshold == "max_sens_spec") %>%
  select(species, model, SORENSEN, TSS, AUC, BOYCE, IMAE) %>%
  pivot_longer(c(SORENSEN, TSS, AUC, BOYCE, IMAE), names_to = "metric") %>%
  mutate(fam = ifelse(metric %in% c("BOYCE", "IMAE"), "ranking (threshold-free)",
                      "binary (max sens/spec)"),
         fam = factor(fam, c("binary (max sens/spec)", "ranking (threshold-free)")))

p_eval <- ggplot(ev, aes(x = model, y = value, color = species, shape = species)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  geom_line(aes(group = species),
            position = position_dodge(width = 0.6),
            linewidth = 0.25, alpha = 0.5) +
  scale_shape(name = "") + scale_color_jco(name = "") +
  scale_y_continuous(labels = scales::label_number(0.01),
                     breaks = c(0.25, 0.5, 0.75, 1)) +
  facet_nested(fam ~ metric, scales = "free_x",
               nest_line = element_line(color = "black", linewidth = .2)) +
  theme_light(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", size = 8,
                                  margin = margin(0.4, 0.4, 0, 0, unit = "mm")),
        legend.position = "bottom",
        legend.text = element_text(margin = margin(0, 0, 0, -1.5, unit = "mm")),
        legend.margin = margin(0, 0, 0, 0, unit = "mm"),
        legend.box.margin = margin(-3.75, 0, -1.5, -4, unit = "mm"),
        legend.box = "vertical",
        panel.spacing.x = unit(0.5, "mm"), panel.spacing.y = unit(1.35, "mm"),
        plot.margin = margin(0, 0.5, 0, 0, unit = "mm"))
ggsave(paste0(OUT, "/stage2_eval.pdf"), p_eval, width = 180, height = 75, units = "mm")

## ---- fig b: suitability maps (the BOYCE-reversal case) ----------------------
sp <- "Cicerbita_alpina"
S <- readRDS(paste0(POC, "/stage2_", sp, "/suitability_surfaces.rds"))
lf <- st_read("/p/projects/isimip/isimip/dannellq/SDMs/Full/C1_F10/1_Inputs/1_Occurrences/subset_repro.gpkg",
              quiet = TRUE) |> st_transform(4326)
pts <- lf[lf$Art_wiss == "Cicerbita alpina", ]
pts <- st_transform(pts, terra::crs(S$sa))

stack <- c(unwrap(S$sa), unwrap(S$sb), unwrap(S$c3w), unwrap(S$s4))   # saved via terra::wrap
names(stack) <- c("A (local HTR)", "B (whole-range WC)",
                  "C3w (meanw N-SDM)", "C4 (B covariate)")
p_map <- ggplot() +
  geom_spatraster(data = stack, maxcell = 2e6) +
  geom_sf(data = pts, size = 0.15, color = "#C42503") +
  scale_fill_viridis_c(name = "", option = "turbo", limits = c(0, 1),
                       na.value = "transparent") +
  facet_wrap(~lyr, ncol = 2) +
  labs(x = NULL, y = NULL) +
  theme_light(base_size = 10) +
  theme(axis.text = element_text(size = 6),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", size = 8),
        legend.position = "bottom",
        legend.margin = margin(-4, 0, 0, 0, unit = "mm"),
        panel.spacing.x = unit(2, "mm"), panel.spacing.y = unit(1, "mm"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, unit = "mm"))
ggsave(paste0(OUT, "/stage2_maps.pdf"), p_map, width = 180, height = 110, units = "mm")
cat("FIGS-DONE ->", OUT, "\n")
