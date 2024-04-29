
sample_pseudoabs_mod <- function (data, x, y, n, method, rlayer, maskval = NULL, calibarea = NULL,
          sp_name = NULL)
{
  . <- ID <- NULL
  if (!any(c("random", "env_const", "geo_const", "geo_env_const",
             "geo_env_km_const") %in% method)) {
    stop("argument 'method' was misused, available methods random, env_const, geo_const, geo_env_const, and geo_env_km_const")
  }
  rlayer <- rlayer[[1]]
  data <- data[, c(x, y)]
  if (!is.null(calibarea)) {
    rlayer <- rlayer %>% terra::crop(., calibarea) %>% terra::mask(.,
                                                                   calibarea)
  }
  if (any(method %in% "random")) {
    cell_samp <- sample_background(data = data, x = x, y = y,
                                   method = "random", n = n, rlayer = rlayer, maskval = maskval)
  }
  if (any(method == "env_const")) {
    if (is.na(method["env"])) {
      stop("Provide a environmental stack/brick variables for env_const method, \ne.g. method = c('env_const', env=somevar)")
    }
    env <- method[["env"]]
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)),
                           as.vector(terra::ext(rlayer)))
      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }
    envp <- inv_bio_mod(e = env, p = data[, c(x, y)])
    envp <- terra::mask(rlayer, envp)
    cell_samp <- sample_background(data = data, x = x, y = y,
                                   method = "random", n = n, rlayer = envp, maskval = maskval)
  }
  if (any(method == "geo_const")) {
    if (!"width" %in% names(method)) {
      stop("Provide a width value for 'geo_const' method, \ne.g. method=c('geo_const', width='50000')")
    }
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    cell_samp <- sample_background(data = data, x = x, y = y,
                                   method = "random", n = n, rlayer = envp, maskval = maskval)
  }
  if (any(method == "geo_env_const")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'geo_env_const' method, \ne.g. method=c('geo_env_const', width='50000', env=somevar)")
    }
    env <- method[["env"]]
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)),
                           as.vector(terra::ext(rlayer)))
      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio_mod(e = env, p = data[, c(x, y)])
    envp <- (envp2 + envp)
    rm(envp2)
    envp <- terra::mask(rlayer, envp)
    cell_samp <- sample_background(data = data, x = x, y = y,
                                   method = "random", n = n, rlayer = envp, maskval = maskval)
  }
  if (any(method == "geo_env_km_const")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'geo_env_km_const' method, \ne.g. method=c('geo_env_const', width='50000', env=somevar)")
    }
    env <- method[["env"]]
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)),
                           as.vector(terra::ext(rlayer)))
      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio_mod(e = env, p = data[, c(x, y)])
    envp <- (envp2 + envp)
    envp <- terra::mask(rlayer, envp)
    rm(envp2)
    if (!is.null(maskval)) {
      if (is.factor(maskval)) {
        maskval <- which(levels(maskval) %in% as.character(maskval))
        rlayer <- rlayer * 1
      }
      filt <- terra::match(rlayer, maskval)
      rlayer <- terra::mask(rlayer, filt)
      rm(filt)
    }
    envp <- terra::mask(rlayer, envp)
    env_changed <- terra::mask(env, envp)
    env_changed <- terra::as.data.frame(env_changed, xy = TRUE)
    env_changed <- stats::na.exclude(env_changed)
    suppressWarnings(km <- stats::kmeans(env_changed, centers = n))
    cell_samp <- km$centers[, 1:2] %>% data.frame()
    val <- terra::extract(envp, cell_samp, method = "simple",
                          xy = TRUE) %>% dplyr::select(-c(ID, x, y))
    cell_samp <- cell_samp %>% dplyr::mutate(val = val[,
                                                       1])
    cell_samp <- cell_samp[!is.na(cell_samp$val), -3]
    cell_samp <- dplyr::tibble(cell_samp)
    cell_samp$pr_ab <- 0
  }
  colnames(cell_samp) <- c(x, y, "pr_ab")
  if (!is.null(sp_name)) {
    cell_samp <- tibble(sp = sp_name, cell_samp)
  }
  return(cell_samp)
}

inv_bio_mod <- function(e, p) {
  if (!methods::is(e, "SpatRaster")) {
    e <- terra::rast(e)
  }
  r <- bio_mod(data = terra::extract(e, p)[-1], env_layer = e)
  r <- (r - terra::minmax(r)[1]) /
    (terra::minmax(r)[2] - terra::minmax(r)[1])
  r <- r <= 0.01 # environmental constrain
  r[which(r[, ] == FALSE)] <- NA
  return(r)
}

bio_mod <- function(data, env_layer) {
  . <- NULL
  if (class(data)[1] != "data.frame") {
    data <- data.frame(data)
  }
  if (!methods::is(env_layer, "SpatRaster")) {
    env_layer <- terra::rast(env_layer)
  }

  data <- na.omit(data)

  result <- env_layer[[1]]
  result[] <- NA

  minv <- apply(data, 2, min)
  maxv <- apply(data, 2, max)
  vnames <- names(data)

  data_2 <- data %>%
    na.omit() %>%
    apply(., 2, sort) %>%
    data.frame()

  if (nrow(data)==1) {data_2 <- t(data_2)}

  rnk <- function(x, y) {
    b <- apply(y, 1, FUN = function(z) sum(x < z))
    t <- apply(y, 1, FUN = function(z) sum(x == z))
    r <- (b + 0.5 * t) / length(x)
    i <- which(r > 0.5)
    r[i] <- 1 - r[i]
    r * 2
  }

  var_df <- terra::as.data.frame(env_layer)
  var_df <- na.omit(var_df)

  k <- (apply(t(var_df) >= minv, 2, all) &
          apply(t(var_df) <= maxv, 2, all))

  for (j in vnames) {
    var_df[k, j] <- rnk(
      data_2[, j],
      var_df[k, j, drop = FALSE]
    )
  }
  var_df[!k, ] <- 0
  res <- apply(var_df, 1, min)
  result[as.numeric(names(res))] <- res
  return(result)
}

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

spThinAlgMod <- function (rec.df.orig, thin.par, reps) 
{
  reduced.rec.dfs <- vector("list", reps)
  DistMat.save <- dist(rec.df.orig, diag = TRUE, upper = TRUE) %>% as.matrix() < thin.par
  diag(DistMat.save) <- FALSE
  DistMat.save[is.na(DistMat.save)] <- FALSE
  SumVec.save <- rowSums(DistMat.save)
  df.keep.save <- rep(TRUE, length(SumVec.save))
  for (Rep in seq_len(reps)) {
    DistMat <- DistMat.save
    SumVec <- SumVec.save
    df.keep <- df.keep.save
    while (any(DistMat) && sum(df.keep) > 1) {
      RemoveRec <- which(SumVec == max(SumVec))
      if (length(RemoveRec) > 1) {
        RemoveRec <- sample(RemoveRec, 1)
      }
      SumVec <- SumVec - DistMat[, RemoveRec]
      SumVec[RemoveRec] <- 0L
      DistMat[RemoveRec, ] <- FALSE
      DistMat[, RemoveRec] <- FALSE
      df.keep[RemoveRec] <- FALSE
    }
    rec.df <- rec.df.orig[df.keep, , drop = FALSE]
    # colnames(rec.df) <- c("x", "y", "Year")
    reduced.rec.dfs[[Rep]] <- rec.df
  }
  reduced.rec.order <- unlist(lapply(reduced.rec.dfs, nrow))
  reduced.rec.order <- order(reduced.rec.order, decreasing = TRUE)
  reduced.rec.dfs <- reduced.rec.dfs[reduced.rec.order]
  return(reduced.rec.dfs)
}
