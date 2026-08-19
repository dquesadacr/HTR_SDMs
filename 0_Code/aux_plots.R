

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
  panel.spacing.x = unit(2, "mm"),
  panel.spacing.y = unit(1, "mm"),
  # panel.grid.minor = element_blank(),
  plot.margin = margin(0.25, 0.25, 0.25, 0.25, unit = "mm"))

prettyZero <- function(l){
    # calculate max.decimals if there are decimal points in the labels, otherwise set to 0
    max.decimals <- ifelse(any(grepl("\\.", l)), 
           max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1,
           0)
    # max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
    lnew = formatC(l, replace.zero = T, zero.print = "0",
        digits = max.decimals, format = "f", preserve.width=T)
    
    # Replace the leading hyphen '-' with the Unicode minus '−' (\u2212)
    # The pattern "^-" ensures we only replace the minus sign at the start of the string
    lnew <- gsub("^-", "\u2010", lnew)

    return(lnew)
}

manual_colors <- c(
                  "#7d5690", 
                  "#6acc68", 
                  "#e39547", 
                  "#34acf8", 
                  "#C42503",
                  "#758796", 
                  "#3f822e", 
                  "#1929C8", 
                  "#ef9fef", 
                  "#C1C1C1", 
                  "#8F0000"
                  )


PlotBRMS2 <- function(brmsobj, tt="Thin", xlabel=TRUE, intercept_text=""){
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
  Block_vec <- c("Predictors", tt, tt, "Approach", "Approach", tt, "Subset", tt)
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
  Post_df$name[which(Post_df$name == "Intercept")] <- intercept_text
  
  Post_df$BLOCK <- factor(Post_df$BLOCK, 
                          levels=c("Intercept","Subset", "Predictors","Approach", 
                                   tt))
  v1 <- eval(parse(text=paste0("c('", Post_df$name %>% unique %>% paste(collapse = "', '"), "')")))
  v1fix <- rev(v1)

  if(all(Post_df$name %>% unique == v1)){
    Post_df <- Post_df %>% 
      mutate(name = fct_relevel(name, v1fix))
  }

  ggplot(Post_df, aes(y = name, x = value*100)) + #
    stat_halfeye(aes(col = Significant), point_size = 1.35, scale=0.925, height=1, slab_fill="#9a9a9a", slab_linewidth=0.35, slab_alpha=0.85,  slab_color="#9a9a9a") + # height=1.25, , linewidth=1, slab_linewidth=1.5 #  slab_color="#9a9a9a",
    facet_nested_wrap(BLOCK~., scales = "free")+
    geom_vline(data = filter(Post_df, BLOCK != "Intercept"), aes(xintercept = 0)) +
    # scale_color_viridis_d(begin = viri_scale[1], end = viri_scale[2]) +
    scale_color_manual(values = c("FALSE" = "#7d5690", "TRUE" = "#6acc68")) +
    guides(color = "none") + 
    theme_light() + 
    labs(x = ifelse(xlabel, "Coefficient Estimate", ""), y = "Coefficients") + 
    gg_theme + 
    theme(plot.margin = margin(0.25,0.25,0.25,1, unit = "mm"),
          axis.title.y = element_blank(),
    )
}
