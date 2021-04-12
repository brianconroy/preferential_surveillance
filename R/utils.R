

equalize_scales <- function(r1, r2){
  
  v1 <- r1[][!is.na(r1[])]
  v2 <- r2[][!is.na(r2[])]
  r_max <- max(v1, v2)
  r_min <- min(v1, v2)
  if (sum(v1 == r_max) == 0){
    v1[length(v1)] <- r_max
  } else{
    v2[length(v1)] <- r_max
  }
  if (sum(v1 == r_min) == 0){
    v1[1] <- r_min
  } else{
    v2[1] <- r_min
  }
  r1[][!is.na(r1[])] <- v1
  r2[][!is.na(r2[])] <- v2
  return(list(r1, r2))
  
}


load_x_ca <- function(factor=NULL){
  
  caPr <- load_prism_pcs()
  if (!is.null(factor)){
    caPr <- aggregate(caPr, fact=factor)
    x_1 <- caPr[[1]][]
    x_1 <- x_1[][!is.na(x_1[])]
    x_2 <- caPr[[2]][]
    x_2 <- x_2[][!is.na(x_2[])]
  } else {
    pc1 <- mask(caPr[[1]], caPr[[2]])
    pc2 <- caPr[[2]]
    x_1 <- pc1[]
    x_1 <- x_1[][!is.na(x_1[])]
    x_2 <- pc2[]
    x_2 <- x_2[][!is.na(x_2[])]
  }
  
  mu_1 <- mean(x_1)
  sd_1 <- sd(x_1)
  mu_2 <- mean(x_2)
  sd_2 <- sd(x_2)
  
  x_1_std <- (x_1 - mu_1)/sd_1
  x_2_std <- (x_2 - mu_2)/sd_2
  x_std <- array(1, c(length(x_1), 1))
  x_std <- cbind(x_std, x_1_std, x_2_std)
  return(x_std)
  
}


assemble_data <- function(rodents, caPr.disc){
  
  coords_all <- cbind(matrix(rodents$Lon_Add_Fix), rodents$Lat_Add_Fix)
  loc.disc <- caPr.disc[[1]]
  cells_all <- cellFromXY(loc.disc, coords_all)
  cells_all <- cells_all[!is.na(cells_all[])]
  cells_obs <- sort(unique(cells_all))
  
  # positive counts at each cell
  rodents_pos <- rodents[rodents$Res == 'POS',]
  coords_pos <- cbind(matrix(rodents_pos$Lon_Add_Fix), rodents_pos$Lat_Add_Fix)
  cells_pos <- cellFromXY(loc.disc, coords_pos)
  counts_pos <- data.frame(table(cells_pos))
  names(counts_pos) <- c('cell', 'count_pos')
  
  # negative counts at each cell
  rodents_neg <- rodents[rodents$Res == 'NEG',]
  coords_neg <- cbind(matrix(rodents_neg$Lon_Add_Fix), rodents_neg$Lat_Add_Fix)
  cells_neg <- cellFromXY(loc.disc, coords_neg)
  counts_neg <- data.frame(table(cells_neg))
  names(counts_neg) <- c('cell', 'count_neg')
  
  # combine counts
  counts_all <- merge(counts_pos, counts_neg, by='cell', all=T)
  counts_all$cell <- as.numeric(as.character(counts_all$cell))
  if (sum(is.na(counts_all$count_pos))){
    counts_all[is.na(counts_all$count_pos),]$count_pos <- 0
  } 
  if (sum(is.na(counts_all$count_neg))){
    counts_all[is.na(counts_all$count_neg),]$count_neg <- 0
  }
  counts_all <- counts_all[with(counts_all, order(cell)),]
  
  # location data
  all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
  loc_raster <- loc.disc
  loc_raster[][!is.na(loc.disc[])] <- 0
  loc_raster[][cells_obs] <- 1
  
  # location covariates
  cov.disc <- caPr.disc
  x1 <- cov.disc[[1]][][all_ids]
  x2 <- cov.disc[[2]][][all_ids]
  x1.standardised <- (x1 - mean(x1))/sd(x1)
  x2.standardised <- (x2 - mean(x2))/sd(x2)
  x <- cbind(1, x1, x2)
  x.standardised <- cbind(1, x1.standardised, x2.standardised)
  locs <- list(
    cells=cells_obs,
    status=1 * c(all_ids %in% cells_obs),  
    coords=xyFromCell(loc.disc, cells_obs),
    raster=loc_raster,
    x.scaled=x.standardised
  )
  locs$ids <- c(1:length(all_ids))[as.logical(locs$status)]
  
  # case data
  cov.disc <- caPr.disc
  x1 <- cov.disc[[1]][][locs$cells]
  x2 <- cov.disc[[2]][][locs$cells]
  x1.standardised <- (x1 - mean(x1))/sd(x1)
  x2.standardised <- (x2 - mean(x2))/sd(x2)
  x <- cbind(1, x1, x2)
  x.standardised <- cbind(1, x1.standardised, x2.standardised)
  
  case.data <- list(
    y=counts_all$count_pos,
    x.standardised=x.standardised,
    x=x,
    p=3
  )
  
  # control data
  ctrl.data <- list(
    y=counts_all$count_neg,
    x.standardised=x.standardised,
    x=x,
    p=3
  )
  
  data <- list(loc=locs, case.data=case.data, ctrl.data=ctrl.data)
  
  return(data)
  
}


reformat_saved_data <- function(data){
  new_locs <- list()
  new_case.data <- list()
  new_ctrl.data <- list()
  for (t in 1:length(data$locs$ids)){
    locs_t <- list(
      status=data$locs$status[[t]],
      cells=data$locs$cells[[t]],
      coords=data$locs$coords[[t]],
      ids=data$locs$ids[[t]]
    )
    new_locs[[t]] <- locs_t
    
    case_t <- list(
      y=data$case.data$y[[t]],
      x.standardised=data$case.data$x.standardised[[t]],
      x=data$case.data$x[[t]],
      p=data$case.data$p[[t]]
    )
    new_case.data[[t]] <- case_t
    
    ctrl_t <- list(
      y=data$ctrl.data$y[[t]],
      x.standardised=data$ctrl.data$x.standardised[[t]],
      x=data$ctrl.data$x[[t]],
      p=data$ctrl.data$p[[t]]
    )
    new_ctrl.data[[t]] <- ctrl_t
  }
  
  new_data <- list(
    locs=new_locs,
    case.data=new_case.data,
    ctrl.data=new_ctrl.data
  )
  return(new_data)
}


add_floor <- function(r){
  
  rvals <- r[][!is.na(r[])]
  rvals[length(rvals)-1] <- 0
  r[][!is.na(r[])] <- rvals
  return(r)
  
}


write_latex_table <- function(df, fname, 
         path="/Users/brianconroy/Documents/research/dataInt/output/"){
  
  df[,] <- lapply(df[, ], as.character)
  tab <- "\\begin{table}[h]
  \\begin{center}
  \\begin{tabular}{l*{5}{c}r}\n"
  header <- paste(names(df), collapse=" & ")
  tab <- paste(tab, header, "\\\\ \n \\hline \n")
  for (i in 1:nrow(df)){
    row_ <- paste(df[i,], collapse=" & ")
    row_ <- paste(row_, '\\\\\n')
    tab <- paste(tab, row_, sep=" ")
  }
  end <- "\\hline
  \\end{tabular}
  \\caption[nrgk]
  {
  nrgk
  }
  \\label{tab:nrgk}
  \\end{center}
  \\end{table}"
  tab <- paste(tab, end, sep=" ")
  path <- paste(path, fname, sep="")
  write(tab, path)
  
}


calc_posterior_lodds <- function(output, x){
  
  n.samp <- nrow(output$samples.beta.ca)
  n.cell <- ncol(output$samples.w)
  lodds_samp <- array(NA, c(n.samp, n.cell))
  
  for (i in 1:n.samp){
    
    beta_ca <- output$samples.beta.ca[i,]
    beta_co <- output$samples.beta.co[i,]
    alpha_ca <- output$samples.alpha.ca[i]
    alpha_co <- output$samples.alpha.co[i]
    w <- output$samples.w[i,]
    
    lodds.i <- x %*% beta_ca + alpha_ca * w - x %*% beta_co - alpha_co * w
    lodds_samp[i,] <- t(lodds.i)
    
  }
  
  return(lodds_samp)
  
}


replace_vals <- function(df, column, val, replacement){
  for (i in 1:nrow(df)){
    if (df[i,column] == val){
      df[i,column] <- replacement
    }
  }
  return(df)
}


ig_var <- function(a, b){
  
  return(b^2/((a-1)^2 * (a-2)))
  
}


g_var <- function(a, b){
  
  return(a*(b^2))

}


n_values <- function(r){
  
  return(length(r[][!is.na(r[])]))
  
}


load_prism_pcs <- function(){
  
  return(raster::stack("data/prism_pcas_ca.grd"))
  
}


#' combine_w
#' 
#' Combines estimated and kriged random effects 
#' according to cell observation status.
#' 
#' @param w_est estimated random effects of observed sites
#' @param w_krige kriged random effects
#' @param status vector of TRUE/FALSE values indicating observation 
#'
#' @return vector of combined estimated and kriged random effects
#' @export
combine_w <- function(w_est, w_krige, status){
  
  w_comb <- c()
  i1 <- 1
  i2 <- 1
  for (i in 1:length(status)){
    if (status[i]){
      w_comb <- c(w_comb, w_est[i1])
      i1 <- i1 + 1
    } else {
      w_comb <- c(w_comb, w_krige[i2])
      i2 <- i2 + 1
    }
  }
  return(w_comb)
  
}


#' view_tr_w
#' 
#' multipanel (4x4) view of random effect traceplots
#'
#' @param samples matrix of mcmc samples
#' @param w_true optional vector of true random effect values
#' @param page view offset
#'
#' @return NULL. Creates plots.
#' @export
view_tr_w <- function(samples, w_true=NULL, page=1){
  
  par(mfrow=c(4, 4))
  for(i in 1:16){
    plot(samples[,page*i], type='l', ylab='w')
    if (!is.null(w_true)){
      abline(h=w_true[page*i], col='2')
    }
  }
  par(mfrow=c(1, 1))
  
}


## arguments
#   output: list, output of a carLeroux function
#   params: list of parameter names to values
summarize <- function(output, params, dic=TRUE){
  
  summary <- list()
  
  counter <- 1
  for (n in names(params)){
    n.summary <- list()
    n.summary[['parameter']] <- n
    n.summary[['true']] <- params[[n]]
    n.summary[['posterior.mean']] <- postMean(n, output)
    n.summary[['posterior.median']] <- postMedian(n, output)
    n.summary[['posterior.sd']] <- postSd(n, output)
    bias.p <- bias(output, n, params[[n]])
    bias.perc.p <- bias(output, n, params[[n]], type='perc')
    n.summary[['bias']] <- bias.p
    n.summary[['percbias']] <- bias.perc.p
    if (dic){
      n.summary[['DIC']] <- output$fit$DIC
      n.summary[['percent deviance']] <- output$fit$percent_dev
    }
    summary[[counter]] <- n.summary
    counter <- counter + 1
  }
  
  if ('rho' %in% names(params)){
    summary[['accept.beta']] <- output$accept[1]
    summary[['accept.phi']] <- output$accept[2]
    summary[['accept.rho']] <- output$accept[3]
  }
  
  return(summary)

}


getSamples <- function(param, output){
  
  if (param == 'beta0'){
    samples <- output$samples.beta[,1]
  } else if (param == 'beta1'){
    samples <- output$samples.beta[,2]
  } else if (param == 'beta2'){
    samples <- output$samples.beta[,3]
  } else if (param == 'beta3'){
    samples <- output$samples.beta[,4]
  } else if (param == 'beta4'){
    samples <- output$samples.beta[,5]
  } else if (param == 'tau2'){
    samples <- output$samples.tau2
  } else if (param == 'tau.sq.1'){
    samples <- output$samples.tau2[,1]
  } else if (param == 'tau.sq.2'){
    samples <- output$samples.tau2[,2]
  } else if (param == 'sigma2'){
    samples <- output$samples.sigma2
  } else if (param == 'rho'){
    samples <- output$samples.rho
  } else if (param == 'rho.1'){
    samples <- output$samples.rho[1,,]
  } else if (param == 'rho.2'){
    samples <- output$samples.rho[2,,]
  } else if (param == 'beta0.loc'){
    samples <- output$samples.beta.l[,1]
  } else if (param == 'beta1.loc'){
    samples <- output$samples.beta.l[,2]
  } else if (param == 'beta2.loc'){
    samples <- output$samples.beta.l[,3]
  } else if (param == 'beta0.cond'){
    samples <- output$samples.beta.c[,1]
  } else if (param == 'beta1.cond'){
    samples <- output$samples.beta.c[,2]
  } else if (param == 'beta2.cond'){
    samples <- output$samples.beta.c[,3]
  } else if (param == 'alpha'){
    samples <- output$samples.alpha
  } else if (param == 'alpha.case'){
    samples <- output$samples.alpha.case
  } else if (param == 'alpha.control'){
    samples <- output$samples.alpha.control
  } else if (param == 'beta0.case'){
    samples <- output$samples.beta.case[,1]
  } else if (param == 'beta1.case'){
    samples <- output$samples.beta.case[,2]
  } else if (param == 'beta2.case'){
    samples <- output$samples.beta.case[,3]
  } else if (param == 'beta0.control'){
    samples <- output$samples.beta.control[,1]
  } else if (param == 'beta1.control'){
    samples <- output$samples.beta.control[,2]
  } else if (param == 'beta2.control'){
    samples <- output$samples.beta.case[,2]
  } 
  
  return(samples)
  
}


postSd <- function(param, output){
  
  samples <- getSamples(param, output)
  return(round(sd(samples), 3))
  
}


postMean <- function(param, output){
  
  samples <- getSamples(param, output)
  return(round(mean(samples), 3))
  
}


postMedian <- function(param, output){
  
  samples <- getSamples(param, output)
  return(round(median(samples), 3)) 
  
}


## arguments
#   output: output of a CAR model
#   param: name of the parameter
#   trueval: true parameter value
#   type: 'abs' or 'perc' (percent)
#   round: rounds bias to 3 digits
bias <- function(output, param, trueval, type='abs', round=TRUE){
  
  samples <- getSamples(param, output)
  bs <- mean(samples) - trueval
  
  if ((type) == 'perc'){
    bs <- 100*bs/trueval
  }
  
  if (round){
    bs <- round(bs, 3)
  }
  
  return(bs)
  
}


# creates traceplots from a list of results,
# i.e. outputs of either carBYM / carLeroux. 
# param specifies which parameter to plot
viewTraces <- function(results, param, line=NULL, titletype='ini'){
  
  for (r in results){
    
    samples <- getSamples(param, r)
    
    if (titletype=='ini'){
      title <- paste('beta0: ', round(r$beta.initial[1], 2),
                    'beta1:', round(r$beta.initial[2], 2),
                    'tau2:', round(r$tau2.initial, 2))
    } else {
      title <- param
    }
    
    plot(samples, type='l', main=title)
    if (!is.null(line)){
      abline(h=line, col='2')
    }
    
  }
  
}


# creates histograms from a list of results,
# i.e. outputs of either carBYM / carLeroux. 
# param specifies which parameter to plot
viewHists <- function(results, param, line=NULL){
  
  
  for (r in results){
    
    samples <- getSamples(param, r) 
    hist(samples, main=param)
    if (!is.null(line)){
      abline(v=line, col='2')
    }
    
  }
}


# arguments
  # x: raster stack to crop and mask
  # ext: extent
  # spdf: spatial polygons data frame defining regions to keep
cropmask <- function(x, ext, spdf){
  
  dat_sub <- crop(x, ext)
  return(mask(dat_sub, spdf))
  
}
