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
  "ggh4x"
)
sapply(package_vec, install.load.package)

## NON-CRAN PACKAGES ----
if("flexsdm" %in% rownames(installed.packages()) == FALSE){ # flexsdm check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("sjevelazco/flexsdm")
}
library(flexsdm)

Spec_ID <- c(cmdArg("id"), cmdArg("s"), cmdArg("species"), cmdArg("ID"))
uniq_sp <- read_sf("./1_Inputs/1_Occurrences/subset_repro.gpkg") %>%
  .$Art_wiss %>% unique

Spec_ID <- str_subset(uniq_sp, fixed(Spec_ID, ignore_case = TRUE))

Spec4folders <- str_replace(Spec_ID, " ", "_")

dic2replace = c("H2000" = "H2k", "Run10"="R10")

temp_list <- list()
for (i in list.files(path = paste0("./1_Inputs/1_Occurrences/", Spec4folders)) %>% 
     str_subset(pattern = "^anti|.pdf|#", negate = TRUE)) {
      tmp <- readRDS(paste0("./1_Inputs/1_Occurrences/", Spec4folders, "/", i))[["Train"]] %>%
        mutate(exp=i %>% str_remove_all(paste0(".rds$"))) %>%
        filter(pr_ab == 1)
      if("Year_i" %in% names(tmp)){
        tmp <- tmp %>% rename(all_of(c(Year="Year_i", Year_i = "Year")))
      }
      temp_list[[i %>% str_remove_all(paste0(".rds$"))]] <- tmp %>% mutate(Year = as.numeric(Year))
    }

wc <- temp_list %>%
  bind_rows() %>%
  separate_wider_delim(exp, delim = "_",
                       names = c("App", "Trim", "Thin", "Dist", "Timespan"), 
                       too_few = "align_start", cols_remove = FALSE) %>%
  replace_na(list(Timespan="", Dist="")) %>%
  mutate(App = str_to_sentence(App) %>% str_replace_all(dic2replace),
         App = fct_relevel(App, c("H2k", "Dec", "R10")),
         Thin = str_remove_all(Thin, "Obs"),
         Trim = str_remove_all(Trim, "Obs"),
         Timespan = fct_relevel(Timespan %>% str_replace_all("100Y", "70Y"),
                                c("", "0Y", "1Y", "5Y", "10Y", "70Y")),
         Dist = fct_relevel(Dist, c("", "1km", "2km")))

wc_1 <- wc %>%
  filter(pr_ab == 1) %>%
  group_by(App, Trim, Thin, Dist, Timespan, Year) %>%
  summarise(Freq = n())

wc_2 <- wc %>% 
  filter(pr_ab == 1) %>%
  group_by(App, Trim, Thin, Dist, Timespan) %>% summarise(Freq = n())

maxfreq <- wc_1 %>% ungroup() %>% summarise(ObsX = max(Freq)) %>% pull(ObsX)
limsY <- wc_1 %>% ungroup() %>% summarise(Ym=min(Year),Yx = max(Year)) %>% as.data.frame()
limsY$Ym <- plyr::round_any(limsY$Ym, 10, floor)
limsY$Yx <- 2015

ObsFreq1 <- ggplot(wc_1 %>% filter(Timespan==""),
                  aes(y= Freq ,
                      x= Year, fill=Freq)) +
  geom_col(width=2) +
  scale_fill_viridis_c(name="", option = "turbo",limits = c(NA, maxfreq)) +
  scale_y_continuous() +
  scale_x_continuous(breaks = seq(limsY$Ym,limsY$Yx, 20), 
                     limits = c(limsY$Ym,limsY$Yx),
                     minor_breaks = seq(limsY$Ym,limsY$Yx, 5)) +
  facet_nested(Thin+Dist~Trim+App, scales="free_y",
               nest_line = element_line(color="black", linewidth = .2)) +
  labs(y= "Yearly observations") +
  theme_light(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=7, angle = 45, hjust = 1),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", margin = margin(0.1,0.1,0.2,0.2, unit = "mm"), size=7),
    legend.position = "bottom",
    legend.key.width = unit(7.5, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(0.5, "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"))

ObsFreq2 <- ggplot(wc_1 %>% filter(Timespan!="", Dist!="", Trim!="NoTrim"),
                   aes(y=Freq, x=Year, fill=Freq)) +
  geom_col(width=2) +
  scale_fill_viridis_c(name="", option = "turbo", limits = c(NA, maxfreq)) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = seq(limsY$Ym,limsY$Yx, 20), 
                     limits = c(limsY$Ym,limsY$Yx),
                     minor_breaks = seq(limsY$Ym,limsY$Yx, 5)) +
  facet_nested_wrap(~App+Dist+Timespan, nrow = 4, dir = "h", 
                    labeller =label_wrap_gen(multi_line = FALSE),
                    nest_line = element_line(color="black", linewidth = .2)) +
  labs(y= "Yearly observations") +
  theme_light(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=7, angle = 45, hjust = 1),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", margin = margin(0.1,0.1,0.2,0.2, unit = "mm"), size=7),
    legend.position = "bottom",
    legend.key.width = unit(7.5, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(0.5, "mm"),
    plot.margin = margin(1, 5.5, 0.5, 0.5, unit = "mm"))

wc_sf <- wc %>%
  st_as_sf(coords = c("x", "y"), crs = "ESRI:31494")

CalibArea <- read_sf("./1_Inputs/3_Calibration_area/Kreis_sud_2km.gpkg") # Calibration Area

ObsMap1 <- ggplot(wc_sf %>% filter(Timespan==""),
                 aes(color= Year)) +
  geom_sf(size=0.35, alpha=0.7) +
  geom_text(data = wc_sf %>% filter(Timespan=="") %>%
              group_by(App, Trim, Thin, Dist, Timespan) %>% summarise(Freq = n()),
            aes(label=paste0("n=" ,Freq)), x= 4640000, y = 5570000, color="black", size=2.25) +
  scale_color_viridis_c(name="", option = "turbo",
                        limits = c(limsY$Ym,limsY$Yx),
                        breaks = seq(limsY$Ym,limsY$Yx, 20)) +
  facet_nested(Thin+Dist~Trim+App,
               nest_line = element_line(color="black", linewidth = .2)) +
  geom_sf(data = CalibArea, color = "black", fill=NA, shape="line") +
  coord_sf() +
  scale_x_continuous(breaks = c(12:14)) +
  scale_y_continuous(breaks = seq(50.4,51.2, 0.4)) +
  theme_light(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", margin = margin(0.1,0.1,0.2,0.2, unit = "mm"), size=7),
    legend.position = "bottom",
    legend.key.width = unit(10, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(0.5, "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"))

ObsMap2 <- ggplot(wc_sf %>% filter(Timespan!="", Dist!="", Trim!="NoTrim"),
                 aes(color= Year)) +
  geom_sf(size=0.35, alpha=0.7) +
  geom_text(data = wc_sf %>% filter(Timespan!="",Dist!="", Trim!="NoTrim") %>%
              group_by(App, Trim, Thin, Dist, Timespan) %>%
              summarise(Freq = n()),
            aes(label=paste0("n=" ,Freq)),
            x= 4640000, y = 5570000, color="black", size=2.25) +
  scale_color_viridis_c(name="", option = "turbo",limits = c(limsY$Ym,limsY$Yx),
                        breaks = seq(limsY$Ym,limsY$Yx, 20)) +
  facet_nested_wrap(~App+Dist+Timespan, nrow = 4, dir = "h",
                    labeller =label_wrap_gen(multi_line = FALSE),
                    nest_line = element_line(color="black", linewidth = .2),) +
  geom_sf(data = CalibArea, color = "black", fill=NA, shape="line") +
  coord_sf() +
  scale_x_continuous(breaks = c(12:14)) +
  scale_y_continuous(breaks = seq(50.4,51.2, 0.4)) +
  theme_light(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color= "black", margin = margin(0.1,0.1,0.2,0.2, unit = "mm"), size=7),
    legend.position = "bottom",
    legend.key.width = unit(10, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    panel.spacing = unit(0.5, "mm"),
    plot.margin = margin(1, 5.5, 0.5, 0.5, unit = "mm"))

joint_plot <- ggarrange(ggarrange(plotlist=list(ObsFreq1, ObsFreq2), nrow = 2,
                                  labels =c("(a)", "(c)"), label.x = -0.015,
                                  font.label = list(size=8, face="plain"), 
                                  common.legend = TRUE, legend = "bottom",
                                  heights = c(.725,1)),
                        ggarrange(plotlist=list(ObsMap1, ObsMap2), nrow = 2,
                                  labels =c("(b)", "(d)"), label.x = -0.015,
                                  font.label = list(size=8, face="plain"), 
                                  common.legend = TRUE, legend = "bottom",
                                  heights = c(0.725,1)),
                        ncol=2, widths = c(1, 1.2)) #heights = c(1,.75)

ggsave(plot= joint_plot,
       filename=paste0("./1_Inputs/1_Occurrences/",Spec4folders, "/thinning.pdf"),
       width=180, height=130, units="mm")
