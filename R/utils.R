

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


calc_posterior_risk <- function(output, x, null_alphas=F){
  
  n.samp <- nrow(output$samples.beta.ca)
  n.cell <- ncol(output$samples.w)
  risk_samp <- array(NA, c(n.samp, n.cell))
  
  for (i in 1:n.samp){
    
    beta_ca <- output$samples.beta.ca[i,]
    beta_co <- output$samples.beta.co[i,]
    alpha_ca <- output$samples.alpha.ca[i]
    alpha_co <- output$samples.alpha.co[i]
    w <- output$samples.w[i,]
    
    if (null_alphas){
      alpha_ca <- 0
      alpha_co <- 0
    }
    
    lodds.i <- x %*% beta_ca + alpha_ca * w - x %*% beta_co - alpha_co * w
    risk_samp[i,] <- t(calc_risk(lodds.i))
    
  }
  
  return(risk_samp)
  
}


calc_posterior_risk_ds <- function(output, x, caPr, caPr.disc, null_alphas=F){
  
  n.samp <- nrow(output$samples.beta.ca)
  
  n.cell_ds <- length(caPr[[1]][][!is.na(caPr[[1]][])])
  cell_ds <- c(1:ncell(caPr[[1]]))[!is.na(values(caPr[[1]]))]
  
  n.cell_lr <- length(caPr.disc[[1]][][!is.na(caPr.disc[[1]][])])
  cells_lr <- c(1:ncell(caPr.disc[[1]]))[!is.na(values(caPr.disc[[1]]))]
  
  risk_samp <- array(NA, c(n.samp, n.cell_ds))
  
  for (i in 1:n.samp){
    
    beta_ca <- output$samples.beta.ca[i,]
    beta_co <- output$samples.beta.co[i,]
    alpha_ca <- output$samples.alpha.ca[i]
    alpha_co <- output$samples.alpha.co[i]
    w <- output$samples.w[i,]
    
    xy_ds <- xyFromCell(caPr[[1]], cell_ds)
    mapped_lr <- cellFromXY(caPr.disc[[1]], xy_ds)
    mapped_lr <- data.frame(cbind(cell_ds, mapped_lr))
    names(mapped_lr) <- c('cell_ds', 'cells_lr')
    
    w2id <- data.frame(cbind(c(1:length(w)), cells_lr))
    names(w2id) <- c('w_id', 'cells_lr')
    wid2ds <- merge(w2id, mapped_lr, by='cells_lr')
    wid2ds <- wid2ds[with(wid2ds, order(cell_ds)),]
    w_ds <- w[wid2ds$w_id]
    
    if (null_alphas){
      alpha_ca <- 0
      alpha_co <- 0
    }
    
    lodds.i <- x %*% beta_ca + alpha_ca * w_ds - x %*% beta_co - alpha_co * w_ds
    risk_samp[i,] <- t(calc_risk(lodds.i))
    
  }
  
  return(risk_samp)
  
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


calc_risk <- function(lodds){
  
  return(exp(lodds)/(1 + exp(lodds)))
  
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



#' summarize_param
#'
#' @param name parameter name
#' @param true true parameter value
#' @param estimated estimated parameter value
#'
#' @return list summarizing estimate and bias
#' @export
#'
#' @examples
summarize_param <- function(name, true, estimated){
  
  summ <- list()
  summ$name <- name
  summ$true <- true
  summ$estimated <- round(estimated, 2)
  summ$bias <- round(estimated - true, 2)
  return(summ)
  
}


#' overlay
#'
#' @param vals vector of ordered cell values
#' @param r raster
#'
#' @return new raster with replaced values
#' @export
#'
#' @examples
overlay <- function(vals, r){
  
  r[][!is.na(r[])] <- vals
  return(r)
  
}


#' combine_w
#' 
#' combines estimated and kriged random effects 
#' according to cell observation status 
#' 
#' @param w_est estimated random effects of observed sites
#' @param w_krige kriged random effects
#' @param status vector of TRUE/FALSE values indicating observation 
#'
#' @return vector of combined estimated and kriged random effects
#' @export
#'
#' @examples
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


#' view_tr
#'
#' @param samples vector of mcmc samples 
#' @param trueval optional true parameter value displayed as horizontal line
#' @param title optional plot title
#' @param ylim optional y axis limits
#'
#' @return
#' @export
#'
#' @examples
view_tr <- function(samples, trueval=NULL, title='', ylim=NULL){
  
  if (is.null(ylim)){
    plot(samples, type='l', main=title)
  } else {
    plot(samples, type='l', main=title, ylim=ylim)
  }
  
  if (!is.null(trueval)){
    abline(h=trueval, col='2')
  }
  
}


#' view_tr_w
#' 
#' multipanel (4x4) view of random effect traceplots
#'
#' @param samples matrix of mcmc samples
#' @param w_true optional vector of true random effect values
#' @param page view offset
#'
#' @return
#' @export
#'
#' @examples
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
