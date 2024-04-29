
plot_preds <- function(ind_p, ens_train, preds, period, posteriori = FALSE, expmt_fold){
  if(posteriori) {
    ens_type <- "predict_posteriori_"
    data2plot <- ind_p
  } else {
    ens_type <- "predict_"
    data2plot <- ind_p$meanw
  }
  dir.create(paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/3_plots/",
                    Spec4folders, "/", ens_type,
                    ens_train, "/", preds),
             recursive = TRUE, showWarnings = FALSE)
  
  plot_ind <- ggplot() +
    geom_spatraster(data = data2plot$max_sens_spec) +
    scale_fill_viridis_d(direction = 1, name ="", begin = 0.2, end = 0.8, na.value = NA, breaks=c("TRUE", "FALSE")) +
    new_scale_fill() +
    geom_spatraster(data = data2plot$meanw) +
    scale_fill_whitebox_c(palette = "soft",na.value = NA,
                          name="", direction=-1, limits=c(0,1)) +
    facet_wrap(~lyr) +
    theme_light(base_size = 10) +
    labs(title= period %>% str_remove("Y")) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(color= "black"),
          legend.position = "bottom",
          legend.key.width = unit(7.5, "mm"),
          legend.key.height = unit(3, "mm"),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(plot_ind, filename=paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/3_plots/",
                                   Spec4folders, "/", ens_type, ens_train, "/",
                                   preds, "/", period %>% str_remove("Y"), ".png"),
         width=15, height=8, units="cm")
}

gif_preds <- function(ens_train, preds, posteriori = FALSE, expmt_fold){
  if(posteriori) {ens_type <- "predict_posteriori_"} else {ens_type <- "predict_"}
  
  if(!str_detect(preds,"H2000")){
    system(paste0("convert -delay 20 -loop 0 ./2_Outputs/", expmt_fold, "/Ensemble/meanw/3_plots/",
                  Spec4folders, "/", ens_type, ens_train, "/", preds, "/*png ",
                  "./2_Outputs/", expmt_fold, "/Ensemble/meanw/3_plots/",
                  Spec4folders, "/", ens_type, ens_train, "/", preds, "/anim.gif"))
  }
}

ens_predict <- function(ensemble_list, predictors_list, plot = TRUE,
                        expmt_fold, preds, fileid) {
  pred_all <- sapply(names(ensemble_list), simplify = FALSE, function(x){
    preds2loop <- names(predictors_list)
    per_ens <- sapply(preds2loop, simplify = FALSE, function(y){
      per_period <- sapply(names(predictors_list[[y]]), simplify = FALSE, function(z){
        print(paste0(z, " - ", y, " - ", x))
        pred <- sdm_predict(
          models = ensemble_list[[x]],
          pred = predictors_list[[y]][[z]],
          thr = "max_sens_spec",
          con_thr = FALSE,
          predict_area = NULL
        )
        if(plot) {plot_preds(ind_p = pred, ens_train = x, preds = y, period = z, expmt_fold = expmt_fold)}
        pred <- sapply(pred, simplify = FALSE, function(x) {x <- terra::wrap(x); return(x)})
        return(pred)
      })
      if(plot) {gif_preds(ens_train = x, preds = y, expmt_fold = expmt_fold)}
      return(per_period)
    })
  })
  # Saving the full list with con and bin in con
  dir.create(paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/1_con/", Spec4folders),
             recursive = TRUE, showWarnings = TRUE)
  saveRDS(pred_all, paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/1_con/", 
                           Spec4folders, "/preds_ens_", fileid, ".rds"))
  return(pred_all)
}

unwrap_lists <- function(x){
  if (is.list(x)) {x <- sapply(x, unwrap_lists, simplify = FALSE)
  } else {x <- terra::unwrap(x)}
}

to_pixels <- function(x){
  if (is.list(x)) {x <- sapply(x, to_pixels, simplify = FALSE)
  } else if (!is.null(x)) {
    freq(x, digits = 1, value = 1) %>% filter(layer==2)}
}

summp_pixels <- function(x, p_inc=0.10){
  if (is.list(x)) {x <- sapply(x, summp_pixels, simplify = FALSE)
  } else if (!is.null(x)) {
    as.data.frame(x$meanw) %>% 
      drop_na() %>%
      # mutate(pgroup = paste0("p", plyr::round_any(meanw, 0.05, floor)*100)) %>%
      mutate(pgroup = plyr::round_any(meanw, p_inc, floor)*100) %>%
      group_by(pgroup) %>% summarise(mean=mean(meanw, na.rm = TRUE), 
                                     n=n()) %>%
      mutate(n_t = sum(n), perc = n/n_t * 100)
  }
}

binPixels2vec <- function(x){
  if (is.list(x)) {x <- sapply(x, binPixels2vec, simplify = FALSE)
  } else if (!is.null(x)) {
    as.data.frame(x$max_sens_spec) %>% 
      drop_na()
  }
}

contPixels2vec <- function(x){
  if (is.list(x)) {x <- sapply(x, contPixels2vec, simplify = FALSE)
  } else if (!is.null(x)) {
    as.data.frame(x$meanw) %>%
      drop_na()
  }
}

binPixels2sum <- function(x){
  if (is.list(x)) {x <- sapply(x, binPixels2sum, simplify = FALSE)
  } else if (!is.null(x)) {
    as.data.frame(x$max_sens_spec) %>% 
      drop_na() %>% sum
  }
}

binPixels2topo <- function(x){
  if (is.list(x)) {x <- sapply(x, binPixels2topo, simplify = FALSE)
  } else if (!is.null(x)) {
    tmp <- DEM*x$max_sens_spec
    tmp[tmp<1e-6] <- NA
    mean(values(tmp), na.rm=TRUE)
  }
}

binPixels2topoPix <- function(x){
  if(any(str_detect(names(x), "orig"))){x <- x[!str_detect(names(x), "orig")]}
  if (is.list(x)) {x <- sapply(x, binPixels2topoPix, simplify = FALSE)
  } else if (!is.null(x)) {
    tmp <- DEM*x$max_sens_spec
    tmp[tmp<1e-6] <- NA
    return(list(Elev=mean(values(tmp), na.rm=TRUE), 
           Pixsum = as.data.frame(x$max_sens_spec) %>% drop_na() %>% sum ))
  }
}

niche_rec <- function(x, y, name = NULL){
  if (is.list(x)) {name <- names(x); x <- sapply(x, function(x) niche_rec(x, name), simplify = FALSE)
  } else if (!is.null(x)) {
    as.data.frame(x$max_sens_spec) %>% 
      drop_na()
  }
}

pred_list_process <- function(x){
  pred_list <- x %>% 
    list_flatten(name_spec = "{outer}-{inner}") %>%
    list_flatten(name_spec = "{outer}-{inner}") %>%
    list_flatten(name_spec = "{outer}-{inner}") %>%
    list_flatten(name_spec = "{outer}-{inner}") %>%
    list_flatten(name_spec = "{outer}-{inner}") %>%
    discard(is.null) %>%
    do.call(rbind,.) %>%
    mutate(Pred = rownames(.)) %>%
    as_tibble() %>%
    separate_wider_delim(cols = Pred, delim = "-", 
                         names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric"), 
                         too_few = "align_start") %>%
    mutate(Posteriori = str_detect(TrainEnsemble, "teriori"),
           Posteriori = case_match(Posteriori, TRUE ~ "Post", FALSE ~ "Reg"),
           Posteriori = factor(Posteriori, c("Reg", "Post")),
           TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
           TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
           PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
    separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach")) %>%
    separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
    separate(TrainApproach, sep = "-", c("TrainApproachPeriod", "TrainApproachObsTrim")) %>%
    separate(PredictApproach, sep = "_", c("ExPreds", "ExSpp")) %>%
    separate(ExPreds, sep = "-", c("ExPredictors", "ExPredsType")) %>%
    separate(ExSpp, sep = "-", c("PredObsPeriod", "PredObsTrim")) %>%
    mutate(TrainApproachObsTrim = case_when(
      TrainApproachPeriod == "H2000" ~ tidyr::replace_na(TrainApproachObsTrim, "NoTrim"),
      TRUE ~ tidyr::replace_na(TrainApproachObsTrim, "ObsTrim")),
      TrainApproachObsTrim = paste0("Train",TrainApproachObsTrim),
      TrainPredsType = case_when(
        !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
        TRUE ~ NA_character_),
      PredObsTrim = case_when(
        PredObsPeriod == "H2000" ~ tidyr::replace_na(PredObsTrim, "NoTrim"),
        TRUE ~ tidyr::replace_na(PredObsTrim, "ObsTrim")),
      PredObsTrim = paste0("Pred",PredObsTrim),
      ExPredsType = case_when(
        !is.na(ExPredsType) ~ paste0("Pred",str_to_sentence(ExPredsType)),
        TRUE ~ NA_character_),
      TrainApproachPeriod = factor(
        paste0("Train", str_to_sentence(TrainApproachPeriod)),
        paste0("Train", c("H2000", "Dec", "Run10"))),
      PredObsPeriod = factor(paste0("Pred", str_to_sentence(PredObsPeriod)),
                             paste0("Pred", c("H2000", "Dec", "Run10"))),
      TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
                               c("WC", "BV", "MI", "B&I")),
      ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
                            c("WC", "BV", "MI", "B&I")),
      Year = str_split_i(Year,"\\.", 1))
  return(pred_list)
}

pred_list_process_grcm <- function(x){
  pred_list <- x %>% 
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    discard(is.null) %>%
    do.call(rbind,.) %>%
    mutate(Pred = rownames(.)) %>%
    as_tibble() %>%
    separate_wider_delim(cols = Pred, delim = "+",
                         names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric"), 
                         too_few = "align_start") %>%
    mutate(TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
           TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
           PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
    separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach")) %>%
    separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
    separate(TrainApproach, sep = "-", c("TrainApproachPeriod", "TrainApproachObsTrim")) %>%
    separate(PredictApproach, sep = "_", 
             c("ExPredictors", "PredsGRCM", "PredsPeriod")) %>%
    separate(PredsGRCM, sep = "-", 
             c("PredsGCM", "PredsRCM", "PredsRCMver", "PredsRealisation", "PredsExp", "PredsDetSto")) %>%
    mutate(TrainApproachObsTrim = case_when(
      TrainApproachPeriod == "H2000" ~ tidyr::replace_na(TrainApproachObsTrim, "NoTrim"),
      TRUE ~ tidyr::replace_na(TrainApproachObsTrim, "ObsTrim")),
      TrainApproachObsTrim = paste0("Train",TrainApproachObsTrim),
      TrainPredsType = case_when(
        !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
        TRUE ~ NA_character_),
      PredsGRCM = paste0(PredsGCM, "-", PredsRCM, "-", PredsRCMver, "-", PredsRealisation),
      PredsExp = case_when(str_detect(PredsExp, "hist") ~ str_to_sentence(PredsExp) %>% str_remove_all("orical"),
                           TRUE ~ str_to_upper(PredsExp)),
      TrainApproachPeriod = factor(
        paste0("Train", str_to_sentence(TrainApproachPeriod)),
        paste0("Train", c("H2000", "Dec", "Run10"))),
      PredsPeriod = factor(paste0("Pred", str_to_sentence(PredsPeriod)),
                           paste0("Pred", c("H2000", "P20", "Dec", "Run10"))),
      TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
                               c("WC", "BV", "MI", "B&I")),
      ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
                            c("WC", "BV", "MI", "B&I")),
      Year = str_split_i(Year,"\\.", 1))
  return(pred_list)
}

niche_list_process_old <- function(x){
  pred_list <- x %>% 
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    discard(is.null) %>%
    do.call(rbind,.) %>% 
    as.data.frame() %>%
    mutate(Pred = rownames(.)) %>%
    unnest(cols = c("SchoenerD", "WarrenI", "HellingerDist")) %>%
    separate_wider_delim(cols = Pred, delim = "+", 
                         names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric", "Threshold"), 
                         too_few = "align_start") %>%
    mutate(TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
           TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
           PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
    separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach")) %>%
    separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
    separate(TrainApproach, sep = "-", c("TrainApproachPeriod", "TrainApproachObsTrim")) %>%
    separate(PredictApproach, sep = "_", 
             c("ExPredictors", "PredsGRCM", "PredsPeriod")) %>%
    separate(PredsGRCM, sep = "-", 
             c("PredsGCM", "PredsRCM", "PredsRCMver", "PredsRealisation", "PredsExp", "PredsDetSto")) %>%
    mutate(TrainApproachObsTrim = case_when(
      TrainApproachPeriod == "H2000" ~ tidyr::replace_na(TrainApproachObsTrim, "NoTrim"),
      TRUE ~ tidyr::replace_na(TrainApproachObsTrim, "ObsTrim")),
      TrainApproachObsTrim = paste0("Train",TrainApproachObsTrim),
      TrainPredsType = case_when(
        !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
        TRUE ~ NA_character_),
      PredsGRCM = paste0(PredsGCM, "-", PredsRCM, "-", PredsRCMver, "-", PredsRealisation),
      PredsExp = case_when(str_detect(PredsExp, "hist") ~ str_to_sentence(PredsExp) %>% str_remove_all("orical"),
                           TRUE ~ str_to_upper(PredsExp)),
      TrainApproachPeriod = factor(
        paste0("Train", str_to_sentence(TrainApproachPeriod)),
        paste0("Train", c("H2000", "Dec", "Run10"))),
      PredsPeriod = factor(paste0("Pred", str_to_sentence(PredsPeriod)),
                           paste0("Pred", c("H2000", "P20", "Dec", "Run10"))),
      TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
                               c("WC", "BV", "MI", "B&I")),
      ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
                            c("WC", "BV", "MI", "B&I")),
      Year = str_split_i(Year,"\\.", 1))
  return(pred_list)
}


niche_list_process_old2 <- function(x){
  pred_list <- x %>% 
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    discard(is.null) %>%
    do.call(rbind,.) %>% 
    as.data.frame() %>%
    mutate(Pred = rownames(.)) %>%
    unnest(cols = c("SchoenerD", "WarrenI", "HellingerDist")) %>%
    separate_wider_delim(cols = Pred, delim = "+",
                         names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric", "Threshold"), 
                         too_few = "align_start") %>%
    mutate(TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
           TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
           PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
    separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach", "TrainTrim", 
                                           "TrainThin", "TrainThinDist", "TrainThinTemp")) %>%
    separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
    separate(PredictApproach, sep = "_",
             c("ExPredictors", "PredsGRCM", "PredsPeriod")) %>%
    separate(PredsGRCM, sep = "-", 
             c("PredsGCM", "PredsRCM", "PredsRCMver", "PredsRealisation", "PredsExp", "PredsDetSto")) %>%
    replace_na(list(TrainThinTemp="NA", TrainThinDist="NA")) %>%
    mutate(
      TrainTrim = paste0("Train",TrainTrim),
      TrainPredsType = case_when(
        !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
        TRUE ~ NA_character_),
      PredsGRCM = paste0(PredsGCM, "-", PredsRCM, "-", PredsRCMver, "-", PredsRealisation),
      PredsExp = case_when(str_detect(PredsExp, "hist") ~ str_to_sentence(PredsExp) %>% str_remove_all("orical"),
                           TRUE ~ str_to_upper(PredsExp)),
      TrainApproach = factor(
        paste0("Train", str_to_sentence(TrainApproach)),
        paste0("Train", c("H2000", "Dec", "Run10"))),
      PredsPeriod = factor(paste0("Pred", str_to_sentence(PredsPeriod)),
                           paste0("Pred", c("H2000", "P20", "Dec", "Run10"))),
      TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
                               c("WC", "BV", "MI", "B&I")),
      ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
                            c("WC", "BV", "MI", "B&I")),
      Year = str_split_i(Year,"\\.", 1),
      TrainThinTemp = fct_relevel(TrainThinTemp, c("NA", "0Y", "5Y", "10Y", "100Y")), #"1Y",
      TrainThinDist = fct_relevel(TrainThinDist %>% str_replace_all("km", " km"), 
                                  c("NA", "1 km", "2 km")))
  return(pred_list)
}


niche_list_process <- function(x){
  dictio_replace_2 <- c("_orig" = "", "_reproj" = "-R", "H2000" = "H2k", "run10"="R10", "dec"="Dec")
  # dictio_replace_preds <- c("worldclim" = "WC", "ind" = "MI", "both"="B&I", "biovars"="BV")

  pred_list <- x %>% 
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    list_flatten(name_spec = "{outer}+{inner}") %>%
    discard(is.null) %>%
    do.call(rbind,.) %>% 
    as.data.frame() %>%
    mutate(Pred = rownames(.)) %>%
    unnest(cols = c("SchoenerD", "WarrenI", "HellingerDist")) %>%
    separate_wider_delim(cols = Pred, delim = "+", 
                         names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric", "Threshold"), 
                         too_few = "align_start") %>%
    mutate(TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
           TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
           PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
    separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach", "TrainTrim", 
                                           "TrainThin", "TrainThinDist", "TrainThinTemp")) %>%
    separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
    separate(PredictApproach, sep = "_", 
             c("ExPredictors", "PredsGRCM", "PredsPeriod")) %>%
    separate(PredsGRCM, sep = "-", 
             c("PredsGCM", "PredsRCM", "PredsRCMver", "PredsRealisation", "PredsExp", "PredsDetSto")) %>%
    replace_na(list(TrainThinTemp="N/A", TrainThinDist="N/A")) %>%
    mutate(
      TrainTrim = TrainTrim %>% str_remove_all("Obs"),
      TrainThin = TrainThin %>% str_remove_all("Obs"),
      TrainPredsType = case_when(
        !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
        TRUE ~ NA_character_),
      PredsGRCM = paste0(PredsGCM, "-", PredsRCM, "-", PredsRCMver, "-", PredsRealisation),
      PredsExp = case_when(str_detect(PredsExp, "hist") ~ str_to_sentence(PredsExp) %>% str_remove_all("orical"),
                           TRUE ~ str_to_upper(PredsExp)),
      TrainApproach = factor(
        TrainApproach %>% str_replace_all(dictio_replace_2),
        c("H2k", "Dec", "R10")),
      PredsPeriod = factor(PredsPeriod %>% str_replace_all(dictio_replace_2),
                           c("H2k", "Dec", "R10")),
      TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
                               c("WC", "BV", "MI", "B&I")),
      ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
                            c("WC", "BV", "MI", "B&I")),
      Year = str_split_i(Year,"\\.", 1),
      TrainThinTemp = fct_relevel(TrainThinTemp %>% str_replace_all("100Y", "70Y"), c("N/A", "0Y", "1Y", "5Y", "10Y", "70Y")),
      TrainThinDist = fct_relevel(TrainThinDist %>% str_replace_all("km", " km"), 
                                  c("N/A", "1 km", "2 km")))
  return(pred_list)
}
