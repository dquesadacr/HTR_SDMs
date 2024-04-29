

gg_theme <- theme_light(base_size = 10) + theme(
  axis.title.x=element_text(size = 9),
  axis.title.y=element_text(size = 9),
  strip.background = element_rect(fill = "white"),
  strip.text = element_text(color= "black", #size=8,
                            margin = margin(0.4,0.4,0.4,0.4, unit = "mm")),
  legend.position = "bottom",
  legend.margin = margin(-4, 0, 0, 0, unit = "mm"),
  legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
  legend.spacing.y = unit(2, "mm"),
  panel.spacing.x = unit(1, "mm"),
  panel.spacing.y = unit(1, "mm"),
  plot.margin = margin(0.25, 0.25, 0.25, 0.25, unit = "mm"))

PlotBRMS2 <- function(brmsobj, tt="Thin"){
  Post_df <- brms::as_draws_df(brmsobj)
  Post_df <- Post_df[ , grep(colnames(Post_df), pattern = "b_")]
  Post_df[,-1] <- apply(Post_df[,-1], MARGIN = 2, FUN = function(x){
    inv_logit_scaled(Post_df$b_Intercept + x) - 
      inv_logit_scaled(Post_df$b_Intercept)
  })
  Post_df$b_Intercept <- inv_logit_scaled(Post_df$b_Intercept)
  Post_df <- Post_df %>%
    pivot_longer(names(Post_df))
  Post_df$name <- gsub(Post_df$name, pattern = "b_", replacement = "")
  Post_df$BLOCK <- NA
  names_vec <- c("PREDICTORS", "Train_THIN", "Train_TRIM", "Pred_TW", "Train_TW", "Pred_T", "Subset", "Train_T")
  Block_vec <- c("Preds", tt, tt, "Approach", "Approach", tt, "Subset", tt)
  for(i in 1:length(names_vec)){
    Post_df$BLOCK[grep(Post_df$name, pattern = names_vec[i])] <- Block_vec[i]
    Post_df$name <- gsub(Post_df$name, pattern = names_vec[i], replacement = "")
  }
  Post_df$name <- gsub(Post_df$name, pattern = "WCMR", replacement = "WC-R")
  Post_df <- Post_df %>% 
    mutate(name=str_replace_all(name, c("TM"="Tr", "PM"="Pr")))
           
  Post_df$BLOCK[Post_df$name == "Intercept"] <- ""
  Post_df$BLOCK[Post_df$name == "n_scale"] <- ""
  Post_df <- Post_df %>% filter(name != "n_scale")
  
  CI <- sapply(unique(Post_df$name), 
               FUN = function(x){
                 x <- Post_df$value[Post_df$name == x]
                 quantile(x, c(.025, 0.975)) # use c(.025, 0.975) for 95%
                 # quantile(x, c(.05, 0.95)) # use c(.025, 0.975) for 95%
                 # quantile(x, c(.005, 0.995)) # use c(.025, 0.975) for 95%
               }
  ) 
  CI_sign <- apply(CI, 2, sign)
  CI_sig <- abs(apply(CI_sign, 2, sum))
  CI_sig <- factor(ifelse(CI_sig == 2, TRUE, FALSE))
  if(all(CI_sig==TRUE)){viri_scale <- c(0.75,0)} else {viri_scale <- c(0,0.75)} # Ghetto fix for the messy colors
  Post_df$Significant <- CI_sig[match(Post_df$name, names(CI_sig))]
  Post_df$name <- str_remove_all(Post_df$name, "as.factor")
  Post_df$BLOCK[which(Post_df$name == "Intercept")] <- "Intercept"
  Post_df$name[which(Post_df$name == "Intercept")] <- ""
  
  Post_df$BLOCK <- factor(Post_df$BLOCK, 
                          levels=c("Intercept","Subset", "Preds","Approach", 
                                   tt))
  v1 <- eval(parse(text=paste0("c('", Post_df$name %>% unique %>% paste(collapse = "', '"), "')")))
  v1fix <- rev(v1)

  if(all(Post_df$name %>% unique == v1)){
    Post_df <- Post_df %>% 
      mutate(name = fct_relevel(name, v1fix))
  }

  ggplot(Post_df, aes(y = name, x = value*100)) + #
    stat_halfeye(aes(col = Significant)) +
    facet_nested_wrap(BLOCK~., scales = "free")+
    geom_vline(data = filter(Post_df, BLOCK != "Intercept"), aes(xintercept = 0)) +
    scale_color_viridis_d(begin = viri_scale[1], end = viri_scale[2]) +
    guides(color = "none") + 
    theme_light() + 
    labs(x = "Coefficient Estimate", y = "Coefficients") + 
    gg_theme + 
    theme(plot.margin = margin(0.25,0.25,0.25,1, unit = "mm"),
          axis.title.y = element_blank(),
    )
}

PlotBRMS3 <- function(brmsobj, tt="Thin & Trim"){
  Post_df <- brms::as_draws_df(brmsobj)
  Post_df <- Post_df[ , grep(colnames(Post_df), pattern = "b_")]
  
  Post_df[,-1] <- apply(Post_df[,-1], MARGIN = 2, FUN = function(x){
    inv_logit_scaled(Post_df$b_Intercept + x) - 
      inv_logit_scaled(Post_df$b_Intercept)
  })
  Post_df$b_Intercept <- inv_logit_scaled(Post_df$b_Intercept)
  
  str_detect(names(Post_df), ":")
  
  Post_df <- Post_df %>%
    pivot_longer(names(Post_df))
  Post_df$name <- gsub(Post_df$name, pattern = "b_", replacement = "")
  
  Post_df$BLOCK <- NA
  names_vec <- c("PREDICTORS", "Train_THIN", "Train_TRIM", "Pred_TW", "Train_TW", "Pred_T", "Subset", "Train_T")
  Block_vec <- c("Preds", tt, tt, "Approach", "Approach", tt, "Subset", tt)
  for(i in 1:length(names_vec)){
    Post_df$BLOCK[grep(Post_df$name, pattern = names_vec[i])] <- Block_vec[i]
    Post_df$name <- gsub(Post_df$name, pattern = names_vec[i], replacement = "")
  }
  Post_df$name <- gsub(Post_df$name, pattern = "WCMR", replacement = "WC-R")
  Post_df <- Post_df %>% 
    mutate(name=str_replace_all(name, c("TM"="Tr", "PM"="Pr")))
  
  Post_df$BLOCK[Post_df$name == "Intercept"] <- ""
  Post_df$BLOCK[Post_df$name == "n_scale"] <- ""
  Post_df <- Post_df %>% filter(name != "n_scale")
  
  CI <- sapply(unique(Post_df$name), 
               FUN = function(x){
                 x <- Post_df$value[Post_df$name == x]
                 quantile(x, c(.025, 0.975)) # use c(.025, 0.975) for 95%
                 # quantile(x, c(.05, 0.95)) # use c(.025, 0.975) for 95%
                 # quantile(x, c(.005, 0.995)) # use c(.025, 0.975) for 95%
               }
  ) 
  CI_sign <- apply(CI, 2, sign)
  CI_sig <- abs(apply(CI_sign, 2, sum))
  CI_sig <- factor(ifelse(CI_sig == 2, TRUE, FALSE))
  if(all(CI_sig==TRUE)){viri_scale <- c(0.75,0)} else {viri_scale <- c(0,0.75)} # Ghetto fix for the messy colors
  Post_df$Significant <- CI_sig[match(Post_df$name, names(CI_sig))]
  Post_df$name <- str_remove_all(Post_df$name, "as.factor")
  Post_df$BLOCK[which(Post_df$name == "Intercept")] <- "Intercept"
  Post_df$name[which(Post_df$name == "Intercept")] <- ""
  
  Post_df$BLOCK <- factor(Post_df$BLOCK, 
                          levels=c("Intercept","Subset", "Preds","Approach", 
                                   tt))
  v1 <- c("","Trim", "Thin", "D2", "Thin:D2")
  v1fix <- rev(v1)
  v2 <- c('', 'MI', 'B&I', 'Dec', 'R10', 
          'D2', 'T0', 'T1', 'T5', 'T10', 'T70', 'D2:T0', 
          'D2:T1', 'D2:T5', 'D2:T10', 'D2:T70')
  v2fix <- rev(v2)
  
  if(all(Post_df$name %>% unique == v1)){
    Post_df <- Post_df %>% 
      mutate(name = fct_relevel(name, v1fix))
  }
  if(all(Post_df$name %>% unique == v2)){
    Post_df <- Post_df %>% 
      mutate(name = fct_relevel(name, v2fix))
  }
  
  ggplot(Post_df, aes(y = name, x = value*100)) + #
    stat_halfeye(aes(col = Significant)) +
    # facet_nested(BLOCK~"Coefficients", scales = "free", 
    #              independent = "x") + 
    facet_nested_wrap(BLOCK~., scales = "free")+
    geom_vline(data = filter(Post_df, BLOCK != "Intercept"), aes(xintercept = 0)) +
    scale_color_viridis_d(begin = viri_scale[1], end = viri_scale[2]) +
    guides(color = "none") + 
    theme_light() + 
    labs(x = "Coefficient Estimate", y = "Coefficients") + 
    gg_theme + 
    theme(plot.margin = margin(0.25,0.25,0.25,1, unit = "mm"),
          axis.title.y = element_blank(),
    )
}
