
plot_preds <- function(sp_ID, ind_p, ens_train, preds, period, expmt_fold, posteriori = FALSE){
  if(posteriori) {
    ens_type <- "predict_posteriori_"
    data2plot <- ind_p
  } else {
    ens_type <- "predict_"
    data2plot <- ind_p$meanw
  }
  dir.create(paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/3_plots/",
                    sp_ID, "/", ens_type,
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
                                   sp_ID, "/", ens_type, ens_train, "/",
                                   preds, "/", period %>% str_remove("Y"), ".png"),
         width=15, height=8, units="cm")
}

gif_preds <- function(sp_ID, ens_train, preds, expmt_fold, posteriori = FALSE){
  if(posteriori) {ens_type <- "predict_posteriori_"} else {ens_type <- "predict_"}
  
  if(!str_detect(preds,"H2000")){
    system(paste0("convert -delay 20 -loop 0 ./2_Outputs/", expmt_fold, "/Ensemble/meanw/3_plots/",
                  sp_ID, "/", ens_type, ens_train, "/", preds, "/*png ",
                  "./2_Outputs/", expmt_fold, "/Ensemble/meanw/3_plots/",
                  sp_ID, "/", ens_type, ens_train, "/", preds, "/anim.gif"))
  }
}

ens_predict <- function(sp_ID, ensemble_list, predictors_list, plot = TRUE,
                        expmt_fold, preds, fileid) {
  pred_all <- sapply(names(ensemble_list), simplify = FALSE, function(x){
    per_ens <- sapply(names(predictors_list), simplify = FALSE, function(y){
      per_period <- sapply(names(predictors_list[[y]]), simplify = FALSE, function(z){
        pred <- sdm_predict(
          models = ensemble_list[[x]],
          pred = predictors_list[[y]][[z]],
          thr = "max_sens_spec",
          con_thr = FALSE,
          predict_area = NULL
        )
        if(plot) {plot_preds(sp_ID = sp_ID, ind_p = pred, ens_train = x, preds = y, period = z, expmt_fold = expmt_fold)}
        pred <- sapply(pred, simplify = FALSE, function(x) {x <- terra::wrap(x); return(x)})
        return(pred)
      })
      if(plot) {gif_preds(sp_ID = sp_ID, ens_train = x, preds = y, expmt_fold = expmt_fold)}
      return(per_period)
    })
  })
  # Saving the full list with con and bin in con
  dir.create(paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/1_con/", sp_ID),
             recursive = TRUE, showWarnings = TRUE)
  saveRDS(pred_all, paste0("./2_Outputs/", expmt_fold, "/Ensemble/meanw/1_con/", 
                           sp_ID, "/preds_ens_", fileid, ".rds"))
  return(pred_all)
}

unwrap_lists <- function(x){
  if (is.list(x)) {x <- sapply(x, unwrap_lists, simplify = FALSE)
  } else {x <- terra::unwrap(x)}
}

eval_df <- function(sp_ID, predictions_list, train_val_list, preds_full, posteriori = FALSE, overwrite = ow_bool){
  if(posteriori) {ens_type <- "eval_posteriori_"} else {ens_type <- "eval_regular_"}
  dir.create(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", sp_ID),
             recursive = TRUE, showWarnings = FALSE)
  if(!overwrite & file.exists(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                        sp_ID, "/", ens_type, preds, ".rds"))){
    eval_all <- readRDS(paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/",
                               sp_ID, "/", ens_type, preds, ".rds"))
  } else {
  eval_all <- sapply(names(predictions_list), simplify = FALSE, function(x){
    per_trama <- sapply(names(predictions_list[[x]]), simplify = FALSE, function(w){

    Trama <- FALSE
      if(str_remove_all(w, preds_full) == (str_remove_all(x, preds_full) %>% str_split_i("_", 1)) ){
        tv_names <- str_subset(names(train_val_list),
                                   str_remove_all(x, ifelse(str_detect(preds_full, "|"), preds_full, w)))
      } else {
        tv_names<- names(train_val_list) %>% str_subset("NoThin") %>%
          str_subset("H2000_Obs", negate=TRUE) %>% str_subset(str_remove_all(w, preds_full) )
        Trama <- TRUE
      }

      per_ens <- sapply(tv_names,simplify = FALSE, function(y){
        per_subset <- sapply(names(train_val_list[[y]]), simplify = FALSE, function(v){
          per_period <- sapply(names(predictions_list[[x]][[w]]), simplify = FALSE, function(z){
          if(z == "H2000"){
            ex_spp <- train_val_list[[y]][[v]]
            if(is.null(ex_spp)) {return()}
          } else {
            if(is.null(train_val_list[[y]][[v]])) {return()}
            ex_spp <- train_val_list[[y]][[v]] %>% filter(Year == z %>% str_remove("Y"))
            if(nrow(ex_spp) == 0 | nrow(ex_spp %>% filter(pr_ab == 0)) == 0) {return()}
          }
          predict_suit <- sdm_extract(
            data = ex_spp,
            x = "x",
            y = "y",
            env_layer = predictions_list[[x]][[w]][[z]]$meanw$meanw,
            variables = NULL,
            filter_na = TRUE
          )
          presence_num <- predict_suit %>% filter(pr_ab == 1) %>% select(meanw) %>% .[[1]] %>% as.numeric
          absence_num <- predict_suit %>% filter(pr_ab == 0) %>% select(meanw) %>% .[[1]] %>% as.numeric

          eval_ens <- sdm_eval(
            p = presence_num,
            a = absence_num
          )

          eval_ens$Year <- z %>% str_remove("Y")
#           print(c(x,y,z,w,v))
          return(eval_ens)
        }) %>% bind_rows()
        per_period$Predictors <- w
        per_period$ObsSubset <- y
        per_period$Subset <- v
        if(Trama) {per_period$Subset <- "Full"}
        per_period$Train <- x
        per_period$Posteriori <- posteriori
        return(per_period)
      }) %>% bind_rows()
    }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
  saveRDS(eval_all, paste0("./2_Outputs/0_Model_performance/Ensemble/meanw/", sp_ID, "/", ens_type, preds, ".rds"))
  }
  return(eval_all)
}

to_pixels <- function(x){
  if (is.list(x)) {x <- sapply(x, to_pixels, simplify = FALSE)
  } else if (!is.null(x)) {
    freq(x, digits = 1, value = 1) %>% filter(layer==2)}
}

ensemble_parse <- function(sp_ID, preds) {
  ens_list <- sapply(list.files(paste0("./2_Outputs/1_Current/Ensemble/meanw/0_models/", sp_ID),
                                paste0("^", preds,"*"), full.names = TRUE), simplify = FALSE, readRDS)
  names(ens_list) <- names(ens_list) %>% str_split_i("/", -1) %>% str_remove(".rds")
  return(ens_list)
}

predictors_parse <- function(preds, expmt_fold) {
  preds_list <- sapply(list.files(paste0("./1_Inputs/2_Predictors/", expmt_fold),
                                  paste0("^",preds,"*"), full.names = TRUE) %>%
                         str_subset("colinvar|full|all", negate = TRUE),
                      simplify = FALSE, function(x){
                        pred <- readRDS(x)
                        if (is(pred, "SpatRaster")) {
                          pred <- list(H2000 = terra::unwrap(pred))
                        } else {
                          pred <- lapply(pred, terra::unwrap)
                        }
                        return(pred)
                      })
  names(preds_list) <- names(preds_list) %>% str_split_i("/", -1) %>% str_remove_all(paste0(".rds"))
  return(preds_list)
}

# summp_pixels <- function(x, p_inc=0.10){
#   if (is.list(x)) {x <- sapply(x, summp_pixels, simplify = FALSE)
#   } else if (!is.null(x)) {
#     as.data.frame(x$meanw) %>%
#       drop_na() %>%
#       # mutate(pgroup = paste0("p", plyr::round_any(meanw, 0.05, floor)*100)) %>%
#       mutate(pgroup = plyr::round_any(meanw, p_inc, floor)*100) %>%
#       group_by(pgroup) %>% summarise(mean=mean(meanw, na.rm = TRUE),
#                                      n=n()) %>%
#       mutate(n_t = sum(n), perc = n/n_t * 100)
#   }
# }

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

# pred_list_process <- function(x){
#   pred_list <- x %>%
#     list_flatten(name_spec = "{outer}-{inner}") %>%
#     list_flatten(name_spec = "{outer}-{inner}") %>%
#     list_flatten(name_spec = "{outer}-{inner}") %>%
#     list_flatten(name_spec = "{outer}-{inner}") %>%
#     list_flatten(name_spec = "{outer}-{inner}") %>%
#     discard(is.null) %>%
#     do.call(rbind,.) %>%
#     mutate(Pred = rownames(.)) %>%
#     as_tibble() %>%
#     separate_wider_delim(cols = Pred, delim = "-",
#                          names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric"),
#                          too_few = "align_start") %>%
#     mutate(Posteriori = str_detect(TrainEnsemble, "teriori"),
#            Posteriori = case_match(Posteriori, TRUE ~ "Post", FALSE ~ "Reg"),
#            Posteriori = factor(Posteriori, c("Reg", "Post")),
#            TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
#            TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
#            PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
#     separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach")) %>%
#     separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
#     separate(TrainApproach, sep = "-", c("TrainApproachPeriod", "TrainApproachObsTrim")) %>%
#     separate(PredictApproach, sep = "_", c("ExPreds", "ExSpp")) %>%
#     separate(ExPreds, sep = "-", c("ExPredictors", "ExPredsType")) %>%
#     separate(ExSpp, sep = "-", c("PredObsPeriod", "PredObsTrim")) %>%
#     mutate(TrainApproachObsTrim = case_when(
#       TrainApproachPeriod == "H2000" ~ tidyr::replace_na(TrainApproachObsTrim, "NoTrim"),
#       TRUE ~ tidyr::replace_na(TrainApproachObsTrim, "ObsTrim")),
#       TrainApproachObsTrim = paste0("Train",TrainApproachObsTrim),
#       TrainPredsType = case_when(
#         !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
#         TRUE ~ NA_character_),
#       PredObsTrim = case_when(
#         PredObsPeriod == "H2000" ~ tidyr::replace_na(PredObsTrim, "NoTrim"),
#         TRUE ~ tidyr::replace_na(PredObsTrim, "ObsTrim")),
#       PredObsTrim = paste0("Pred",PredObsTrim),
#       ExPredsType = case_when(
#         !is.na(ExPredsType) ~ paste0("Pred",str_to_sentence(ExPredsType)),
#         TRUE ~ NA_character_),
#       TrainApproachPeriod = factor(
#         paste0("Train", str_to_sentence(TrainApproachPeriod)),
#         paste0("Train", c("H2000", "Dec", "Run10"))),
#       PredObsPeriod = factor(paste0("Pred", str_to_sentence(PredObsPeriod)),
#                              paste0("Pred", c("H2000", "Dec", "Run10"))),
#       TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
#                                c("WC", "BV", "MI", "B&I")),
#       ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
#                             c("WC", "BV", "MI", "B&I")),
#       Year = str_split_i(Year,"\\.", 1))
#   return(pred_list)
# }
#
# pred_list_process_grcm <- function(x){
#   pred_list <- x %>%
#     list_flatten(name_spec = "{outer}+{inner}") %>%
#     list_flatten(name_spec = "{outer}+{inner}") %>%
#     list_flatten(name_spec = "{outer}+{inner}") %>%
#     list_flatten(name_spec = "{outer}+{inner}") %>%
#     list_flatten(name_spec = "{outer}+{inner}") %>%
#     discard(is.null) %>%
#     do.call(rbind,.) %>%
#     mutate(Pred = rownames(.)) %>%
#     as_tibble() %>%
#     separate_wider_delim(cols = Pred, delim = "+",
#                          names = c("TrainEnsemble", "TrainPredictors", "PredictApproach", "Year", "Metric"),
#                          too_few = "align_start") %>%
#     mutate(TrainEnsemble = str_remove_all(TrainEnsemble, "ens_|posteriori_"),
#            TrainPredictors = str_replace_all(TrainPredictors, dictio_replace),
#            PredictApproach = str_replace_all(PredictApproach, dictio_replace)) %>%
#     separate(TrainPredictors, sep = "_", c("TrainPreds", "TrainApproach")) %>%
#     separate(TrainPreds, sep = "-", c("TrainPredictors", "TrainPredsType")) %>%
#     separate(TrainApproach, sep = "-", c("TrainApproachPeriod", "TrainApproachObsTrim")) %>%
#     separate(PredictApproach, sep = "_",
#              c("ExPredictors", "PredsGRCM", "PredsPeriod")) %>%
#     separate(PredsGRCM, sep = "-",
#              c("PredsGCM", "PredsRCM", "PredsRCMver", "PredsRealisation", "PredsExp", "PredsDetSto")) %>%
#     mutate(TrainApproachObsTrim = case_when(
#       TrainApproachPeriod == "H2000" ~ tidyr::replace_na(TrainApproachObsTrim, "NoTrim"),
#       TRUE ~ tidyr::replace_na(TrainApproachObsTrim, "ObsTrim")),
#       TrainApproachObsTrim = paste0("Train",TrainApproachObsTrim),
#       TrainPredsType = case_when(
#         !is.na(TrainPredsType) ~ paste0("Train",str_to_sentence(TrainPredsType)),
#         TRUE ~ NA_character_),
#       PredsGRCM = paste0(PredsGCM, "-", PredsRCM, "-", PredsRCMver, "-", PredsRealisation),
#       PredsExp = case_when(str_detect(PredsExp, "hist") ~ str_to_sentence(PredsExp) %>% str_remove_all("orical"),
#                            TRUE ~ str_to_upper(PredsExp)),
#       TrainApproachPeriod = factor(
#         paste0("Train", str_to_sentence(TrainApproachPeriod)),
#         paste0("Train", c("H2000", "Dec", "Run10"))),
#       PredsPeriod = factor(paste0("Pred", str_to_sentence(PredsPeriod)),
#                            paste0("Pred", c("H2000", "P20", "Dec", "Run10"))),
#       TrainPredictors = factor(TrainPredictors %>% str_replace_all(dictio_replace_preds),
#                                c("WC", "BV", "MI", "B&I")),
#       ExPredictors = factor(ExPredictors %>% str_replace_all(dictio_replace_preds),
#                             c("WC", "BV", "MI", "B&I")),
#       Year = str_split_i(Year,"\\.", 1))
#   return(pred_list)
# }

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
