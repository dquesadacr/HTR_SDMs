
spt_thin <- function(xyt.df, dist.filt, timespan, reps=20) {
  reduced.rec.dfs <- vector("list", reps)

  for (Rep in 1:reps) {
    df <- xyt.df
    iter <- TRUE

    while(iter) {
      # Calculate distance matrix
      # Filter based on spatial distance
      DistMat <- dist(df[, c("x", "y")], diag = TRUE, upper = TRUE) %>%
        as.matrix() < dist.filt

      # Filter based on temporal distance
      YearVec <- df[, "Year"]
      YearDistMat <- abs(outer(YearVec, YearVec, "-")) <= timespan
      DistMat <- DistMat & YearDistMat

      colnames(DistMat) <- rownames(DistMat) <- NULL
      SumVec <- rowSums(DistMat)

      Year.table <- table(YearVec) %>% sort(decreasing = TRUE)
      YT.years <- names(Year.table) %>% as.numeric()

      if(length(which(SumVec == max(SumVec))) > 1) {
        ind2check <- sample(which(SumVec == max(SumVec)), size = 1)
      } else {ind2check <- which(SumVec == max(SumVec))}

      # Logic to decide which to remove
      candidates <- data.frame(Ind=which(DistMat[,ind2check]),
                               Year = YearVec[which(DistMat[,ind2check])]) %>%
        mutate(YT = match(Year, YT.years)) %>%
        group_by(Year) %>%
        mutate(nPerYear = n()) %>%
        arrange(desc(nPerYear),YT) %>%
        ungroup()

      if(any(candidates$nPerYear > 1)){
        df <-  df[-(candidates %>% filter(nPerYear == max(nPerYear)) %>%
                      slice_sample(n = -1) %>% pull(Ind)),]
      } else if (max(SumVec)>1 & all(candidates$nPerYear == 1)) {
        # break
        df <-  df[-(candidates %>% slice_min(YT, n = -1) %>% pull(Ind)),]
      } else {iter <- FALSE}
    }
    rec.df <- df
    reduced.rec.dfs[[Rep]] <- rec.df
  }
  reduced.rec.order <- unlist(lapply(reduced.rec.dfs, nrow))
  reduced.rec.order <- order(reduced.rec.order, decreasing = TRUE)
  reduced.rec.dfs <- reduced.rec.dfs[reduced.rec.order]
  return(reduced.rec.dfs)
}
