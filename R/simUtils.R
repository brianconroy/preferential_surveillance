library(jsonlite)


get_gamma_prior <- function(prior_mean, prior_var){
  
  shape <- prior_mean^2/prior_var
  scale <- prior_var/prior_mean
  return(c(shape, scale))
  
}


get_igamma_prior <- function(prior_mean, prior_var){
  
  scale <- prior_mean^3 * (1/prior_var + 1/prior_mean^2)
  shape <- scale/prior_mean + 1
  return(c(shape, scale))
  
}


summarize_mcmc_pscc <- function(output, model_desc){
  
  rows <- list()
  
  rows[[1]] <- list(
    model=model_desc,
    parameter='n.sample',
    value=output$n.sample
  )
  rows[[2]] <- list(
    model=model_desc,
    parameter='burnin',
    value=output$burnin
  )
  rows[[3]] <- list(
    model=model_desc,
    parameter='proposal.sd.theta',
    value=output$proposal.sd.theta
  )
  rows[[4]] <- list(
    model=model_desc,
    parameter='L (w)',
    value=output$L_w
  )
  rows[[5]] <- list(
    model=model_desc,
    parameter='L (beta case)',
    value=output$L_ca
  )
  rows[[6]] <- list(
    model=model_desc,
    parameter='L (beta control)',
    value=output$L_co
  )
  rows[[7]] <- list(
    model=model_desc,
    parameter='L (alpha case)',
    value=output$L_a_ca
  )
  rows[[8]] <- list(
    model=model_desc,
    parameter='L (alpha control)',
    value=output$L_a_co
  )
  rows[[9]] <- list(
    model=model_desc,
    parameter='m',
    value=1000
  )
  rows[[10]] <- list(
    model=model_desc,
    parameter='target acceptance',
    value=0.65
  )
  
  return(rows)
  
}


calc_log_odds_output <- function(output, true_params){
  
  location_indicators <- true_params$location_indicators
  x_standard <- f(location_indicators)
  
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  
  lodds.ps <- x_standard %*% beta_ca_h + alpha_ca_h * w.hat - x_standard %*% beta_co_h - alpha_co_h * w.hat
  return(lodds.ps)
  
}


#' table_params
#'
#' @param outputs (list) mcmc outputs
#' @param sampling (character) strength of sampling
#' @param prevalence (character) disease prevalence ("low", "medium", "high")
#'
#' @return (list) list of lists containing parameter
#' estimates, true values and biases of different models
#' @export
#'
#' @examples
table_params <- function(outputs, sampling, prevalence){
  
  output_ps <- get_output(outputs, sampling, prevalence, 'prefSampleGpCC')
  output_sp_ca <- get_output(outputs, sampling, prevalence, 'spatial_poisson_case')
  output_sp_co <- get_output(outputs, sampling, prevalence, 'spatial_poisson_ctrl')
  betas <- load_params(paste('estimates_poisson_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
  betas$description <- "Poisson Regression"
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  
  rows <- list()
  shared <- c("Beta 0 (case)", "Beta 1 (case)", "Beta 2 (case)", 
              "Beta 0 (control)", "Beta 1 (control)", "Beta 2 (control)")
  counter <- 1
  for (p in shared){
    rows[[counter]] <- make_row(
      prevalence,
      sampling,
      "PS",
      p,
      true_params,
      output_ps
    )
    rows[[counter + 1]] <- make_row(
      prevalence,
      sampling,
      "SP",
      p,
      true_params,
      output_sp_ca
    )
    rows[[counter + 2]] <- make_row(
      prevalence,
      sampling,
      "PR",
      p,
      true_params,
      betas
    )
    counter <- counter + 4
  }
  
  ps <- c("Alpha (case)", "Alpha (control)", "Theta", "Phi")
  for (p in ps){
    rows[[counter]] <- make_row(
      prevalence,
      sampling,
      "PS",
      p,
      true_params,
      output_ps
    )
    counter <- counter + 1
  }
  
  return(rows)
}


make_row <- function(prevalence, sampling, model, parameter, true_params, output){
  
  est <- get_estimate(output, parameter)
  true <- get_true_val(true_params, parameter)
  row <- list(
    prevalence=prevalence,
    sampling=sampling,
    model=model,
    parameter=parameter,
    estimate=est,
    true=true,
    bias=round(est-true, 3)
  )
  return(row)
  
}


make_row_wide <- function(outputs, sampling, prevalence, param, output_tag, model_name){
  row <- list(
    Sampling=sampling,
    Model=model_name,
    Parameter=param
  )
  for (prevalence in c("low", "medium", "high")){
    if (output_tag != "poisson"){
      output <- get_output(outputs, sampling, prevalence, output_tag)
    } else {
      output <- load_params(paste('estimates_poisson_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
      output$description <- "Poisson Regression"
    }
    est <- get_estimate(output, param)
    true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
    true <- get_true_val(true_params, param)
    row[prevalence] <- paste(est, " (", round(est-true, 3), ")", sep="")
  }
  return(row)
}


table_params_wide <- function(outputs, sampling){
  
  shared <- c("Beta 0 (case)", "Beta 1 (case)", "Beta 2 (case)", 
              "Beta 0 (control)", "Beta 1 (control)", "Beta 2 (control)")
  rows <- list()
  counter <- 1
  for (param in shared){
    rows[[counter]] <- make_row_wide(outputs, sampling, prevalence, param, 'prefSampleGpCC', 'PrefSample')
    if (grepl("case", param)){
      rows[[counter+1]] <- make_row_wide(outputs, sampling, prevalence, param, 'spatial_poisson_case', 'SpatPoisson')
    } else {
      rows[[counter+2]] <- make_row_wide(outputs, sampling, prevalence, param, 'spatial_poisson_ctrl', 'SpatPoisson')
    }
    rows[[counter+3]] <- make_row_wide(outputs, sampling, prevalence, param, 'poisson', 'Poisson')
    counter <- counter + 4
  }
  return(ldply(rows, 'data.frame'))
  
}


make_row_summary <- function(outputs, sampling, prevalence, param, output_tag, model_name){
  row <- list(
    Sampling=sampling,
    Model=model_name,
    Parameter=param
  )
  for (prevalence in c("low", "medium", "high")){
    if (output_tag != "poisson"){
      output <- get_output(outputs, sampling, prevalence, output_tag)
    } else {
      output <- load_params(paste('estimates_poisson_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
      output$description <- "Poisson Regression"
    }
    est <- get_estimate(output, param)
    true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
    true <- get_true_val(true_params, param)
    row[prevalence] <- round(est-true, 3)
  }
  return(row)
}


table_params_summary <- function(outputs, sampling){
  
  shared <- c("Beta 0 (case)", "Beta 1 (case)", "Beta 2 (case)", 
              "Beta 0 (control)", "Beta 1 (control)", "Beta 2 (control)")
  rows <- list()
  counter <- 1
  for (param in shared){
    rows[[counter]] <- make_row_summary(outputs, sampling, prevalence, param, 'prefSampleGpCC', 'PrefSample')
    if (grepl("case", param)){
      rows[[counter+1]] <- make_row_summary(outputs, sampling, prevalence, param, 'spatial_poisson_case', 'SpatPoisson')
    } else {
      rows[[counter+2]] <- make_row_summary(outputs, sampling, prevalence, param, 'spatial_poisson_ctrl', 'SpatPoisson')
    }
    rows[[counter+3]] <- make_row_summary(outputs, sampling, prevalence, param, 'poisson', 'Poisson')
    counter <- counter + 4
  }
  return(ldply(rows, 'data.frame'))
  
}


get_true_val <- function(true_params, parameter){
  
  if (parameter == "Beta 0 (case)"){
    return(true_params$beta.case[1])
  } else if (parameter == "Beta 1 (case)"){
    return(true_params$beta.case[2])
  } else if (parameter == "Beta 2 (case)"){
    return(true_params$beta.case[3])
  } else if (parameter == "Beta 0 (control)"){
    return(true_params$beta.ctrl[1])
  } else if (parameter == "Beta 1 (control)"){
    return(true_params$beta.ctrl[2])
  } else if (parameter == "Beta 2 (control)"){
    return(true_params$beta.ctrl[3])
  } else if (parameter == "Alpha (case)"){
    return(true_params$Alpha.case)
  } else if (parameter == "Alpha (control)"){
    return(true_params$Alpha.ctrl)
  } else if (parameter == "W"){
    return(true_params$W)
  } else if (parameter == "Theta"){
    return(true_params$Theta)
  } else if (parameter == "Phi"){
    return(true_params$Phi)
  }
  
}


get_estimate <- function(output, parameter){
  
  if (grepl("prefSampleGpCC", output$description)){
    type <- "PS"
  } else if (grepl("Poisson Regression", output$description)){
    type <- "PR"
  } else {
    type <- "SP"
  }
  
  if (parameter == "Beta 0 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,1]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,1]
    } else {
      target_samples <- output$case[1]
    }
  } else if (parameter == "Beta 1 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,2]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,2]
    } else {
      target_samples <- output$case[2]
    }
  } else if (parameter == "Beta 2 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,3]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,3]
    } else {
      target_samples <- output$case[3]
    }
  } else if (parameter == "Beta 0 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,1]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,1]
    } else {
      target_samples <- output$ctrl[1]
    }
  } else if (parameter == "Beta 1 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,2]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,2]
    } else {
      target_samples <- output$ctrl[2]
    }
  } else if (parameter == "Beta 2 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,3]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,3]
    } else {
      target_samples <- output$ctrl[3]
    }
  } else if (parameter == "Alpha (case)"){
    target_samples <- output$samples.alpha.ca
  } else if (parameter == "Alpha (control)"){
    target_samples <- output$samples.alpha.co
  } else if (parameter == "Theta"){
    target_samples <- output$samples.theta
  } else if (parameter == "Phi"){
    target_samples <- output$samples.phi
  }
  
  est <- mean(target_samples)
  return(round(est, 3))
  
}


get_pad <- function(val, typ){
  if (typ == 'lower'){
    if (val < 0){
      return(1.15 * val)
    } else if (val == 0) {
      return(-0.15)
    } else {
      return(0.75 * val)
    } 
  } else {
    if (val < 0){
      return(0.75 * val)
    } else if (val == 0) {
      return(0.15)
    } else {
      return(1.15 * val)
    } 
  }
}


padded_plot <- function(y, trueval, ylab='', title=''){
  
  ymax <- max(y)
  ymin <- min(y)
  if (trueval < ymin){
    lb <- get_pad(trueval, 'lower')
    ub <- get_pad(ymax, 'upper')
    plot(y, typ='l', ylab=ylab, ylim=c(lb, ub), main=title); abline(h=trueval, col='2')
  } else if (trueval > ymax){
    ub <- get_pad(trueval, 'upper')
    lb <- get_pad(ymin, 'lower')
    plot(y, typ='l', ylab=ylab, ylim=c(lb, ub), main=title); abline(h=trueval, col='2')
  } else{
    plot(y, typ='l', ylab=ylab, main=title); abline(h=trueval, col='2')
  }
  
}


padded_plot2 <- function(y, ylab=''){
  
  ymax <- max(y)
  ymin <- min(y)
  if (ymax > 0){
    ub <- 2 * ymax
    if (ymin > 0){
      lb <- 0
    } else{
      lb <- 2 * ymin
    }
  } else {
    lb <- 2 * ymin
    ub <- 2 * abs(ub)
  }
  plot(y, ylim=c(lb, ub), type='l')
  
}


#' plot_traces
#' 
#' plots traceplots for the preferential sampling model
#'
#' @param outputs (list) mcmc preferential sampling outputs
#' @param sampling (character) strength of preferential sampling
#' @param prevalence (character) disease prevalence
#'
#' @return
#' @export
#'
#' @examples
plot_traces <- function(outputs, sampling, prevalence){
  
  output <- get_output(outputs, sampling, prevalence, 'prefSampleGpCC')
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  
  par(mfrow=c(3,4))
  padded_plot(output$samples.beta.ca[,1], ylab='Beta 0 (case)', true_params$beta.case[1])
  padded_plot(output$samples.beta.ca[,2], ylab='Beta 1 (case)', true_params$beta.case[2])
  padded_plot(output$samples.beta.ca[,3], ylab='Beta 2 (case)', true_params$beta.case[3])
  padded_plot(output$samples.beta.co[,1], ylab='Beta 0 (control)', true_params$beta.ctrl[1])
  padded_plot(output$samples.beta.co[,2], ylab='Beta 1 (control)', true_params$beta.ctrl[2])
  padded_plot(output$samples.beta.co[,3], ylab='Beta 2 (control)', true_params$beta.ctrl[3])
  padded_plot(output$samples.alpha.ca, ylab='Alpha (case)', true_params$Alpha.case[1])
  padded_plot(output$samples.alpha.co, ylab='Alpha (control)', true_params$Alpha.ctrl[1])
  padded_plot(output$samples.theta, ylab='Range', true_params$Theta[1])
  padded_plot(output$samples.phi, ylab='Marginal Variance', true_params$Phi)
  par(mfrow=c(1,1))
  
}


plot_traces_general <- function(output, true_params=NULL){
  
  par(mfrow=c(3,4))
  if (!is.null(true_params)){
    padded_plot(output$samples.beta.ca[,1], ylab='Beta 0 (case)', true_params$beta.case[1])
    padded_plot(output$samples.beta.ca[,2], ylab='Beta 1 (case)', true_params$beta.case[2])
    padded_plot(output$samples.beta.ca[,3], ylab='Beta 2 (case)', true_params$beta.case[3])
    padded_plot(output$samples.beta.co[,1], ylab='Beta 0 (control)', true_params$beta.ctrl[1])
    padded_plot(output$samples.beta.co[,2], ylab='Beta 1 (control)', true_params$beta.ctrl[2])
    padded_plot(output$samples.beta.co[,3], ylab='Beta 2 (control)', true_params$beta.ctrl[3])
    padded_plot(output$samples.alpha.ca, ylab='Alpha (case)', true_params$Alpha.case[1])
    padded_plot(output$samples.alpha.co, ylab='Alpha (control)', true_params$Alpha.ctrl[1])
    padded_plot(output$samples.theta, ylab='Range', true_params$Theta[1])
    padded_plot(output$samples.phi, ylab='Marginal Variance', true_params$Phi)
  } else {
    plot(output$samples.beta.ca[,1], ylab='Beta 0 (case)', type='l')
    plot(output$samples.beta.ca[,2], ylab='Beta 1 (case)', type='l')
    plot(output$samples.beta.ca[,3], ylab='Beta 2 (case)', type='l')
    plot(output$samples.beta.co[,1], ylab='Beta 0 (control)', type='l')
    plot(output$samples.beta.co[,2], ylab='Beta 1 (control)', type='l')
    plot(output$samples.beta.co[,3], ylab='Beta 2 (control)', type='l')
    plot(output$samples.alpha.ca, ylab='Alpha (case)', type='l')
    plot(output$samples.alpha.co, ylab='Alpha (control)', type='l')
    plot(output$samples.theta, ylab='Range', type='l')
    plot(output$samples.phi, ylab='Marginal Variance', type='l')
  }
  par(mfrow=c(1,1))
  
}


get_output <- function(outputs, sampling, prevalence, type){
  
  output_target <- list()
  tag <- paste(sampling, prevalence, sep="_")
  for (o in outputs){
    if (grepl(tag, o$description) & grepl(type, o$description)){
      output_target <- o
      break
    }
  }
  return(output_target)
  
}


get_output_general <- function(outputs, tag){
  
  output_target <- list()
  for (o in outputs){
    if (grepl(tag, o$description)){
      output_target <- o
      break
    }
  }
  return(output_target)
  
}


calc_log_odds <- function(outputs, sampling, prevalence){
  
  output_target <- list()
  tag <- paste(sampling, prevalence, sep="_")
  for (o in outputs){
    if (grepl(tag, o$description) & grepl('prefSampleGpCC', o$description)){
      output_target <- o
      break
    }
  }
  
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  
  w.hat <- colMeans(output_target$samples.w)
  beta_ca_h <- colMeans(output_target$samples.beta.ca)
  beta_co_h <- colMeans(output_target$samples.beta.co)
  alpha_ca_h <- mean(output_target$samples.alpha.ca)
  alpha_co_h <- mean(output_target$samples.alpha.co)
  
  lodds.ps <- x_standard %*% beta_ca_h + alpha_ca_h * w.hat - x_standard %*% beta_co_h - alpha_co_h * w.hat
  return(lodds.ps)
  
}


calc_log_odds_true <- function(sampling, prevalence){
  
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  beta.case <- true_params$beta.case
  beta.ctrl <- true_params$beta.ctrl
  Alpha.case <- true_params$Alpha.case
  Alpha.ctrl <- true_params$Alpha.ctrl
  W <- true_params$W
  
  lodds.true <- x_standard %*% beta.case + Alpha.case * W - x_standard %*% beta.ctrl - Alpha.ctrl * W
  return(lodds.true)
  
}


calc_log_odds_true_general <- function(true_params){
  
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  beta.case <- true_params$beta.case
  beta.ctrl <- true_params$beta.ctrl
  Alpha.case <- true_params$Alpha.case
  Alpha.ctrl <- true_params$Alpha.ctrl
  W <- true_params$W
  
  lodds.true <- x_standard %*% beta.case + Alpha.case * W - x_standard %*% beta.ctrl - Alpha.ctrl * W
  return(lodds.true)
  
}


#' calc_log_odds_sp
#' 
#' calculates the log odds from the spatial poisson models
#'
#' @param outputs 
#' @param sampling 
#' @param prevalence 
#'
#' @return
#' @export
#'
#' @examples
calc_log_odds_sp <- function(outputs, sampling, prevalence){
  
  tag <- paste(sampling, prevalence, sep='_')
  output_ca <- list()
  output_co <- list()
  for (o in outputs){
    if (grepl(tag, o$description) & grepl('spatial_poisson_case', o$description)){
      output_ca <- o
    } else if (grepl(tag, o$description) & grepl('spatial_poisson_ctrl', o$description)){
      output_co <- o
    }
  }
  
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  location_ids <- true_params$location_ids
  x_standard <- load_x_standard(location_indicators)
  
  w.hat_spca <- colMeans(output_ca$samples.w)
  beta_ca_sp <- colMeans(output_ca$samples.beta)
  kriged_w_ca <- load_output(paste('output.krige_ca_prefSampleGpCC_', sampling, '_',prevalence, '.json', sep=''))
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, location_indicators)
  
  w.hat_spco <- colMeans(output_co$samples.w)
  beta_co_sp <- colMeans(output_co$samples.beta)
  kriged_w_co <- load_output(paste('output.krige_co_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, location_indicators)
  
  lodds <- x_standard %*% beta_ca_sp + w_ca_est - x_standard %*% beta_co_sp - w_co_est
  return(lodds)
  
}


#' calc_log_odds_pr
#' 
#' calculates the log odds from the poisson regression models
#'
#' @param sampling 
#' @param prevalence 
#'
#' @return
#' @export
#'
#' @examples
calc_log_odds_pr <- function(sampling, prevalence){
  
  betas <- load_params(paste('estimates_poisson_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
  beta_ca_r <- betas$case
  beta_co_r <- betas$ctrl
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  lodds <- x_standard %*% beta_ca_r - x_standard %*% beta_co_r
  return(lodds)
  
}


#' save_output
#'
#' @param output (list) output of an MCMC function
#' @param fname (character) file name
#'
#' @return writes output as json file
#'
save_output <- function(output, fname, dst="/Users/brianconroy/Documents/research/dataInt/output/"){
  
  path <- paste(dst, fname, sep="")
  write(toJSON(output), path)
  
}


#' load_output
#'
#' @param file (character) name of output JSON file to be loaded
#'
#' @return (list) loaded output
#' @export
#'
#' @examples load_output("prefSampleGpCC_med.json")
load_output <- function(fname, src="/Users/brianconroy/Documents/research/dataInt/output/"){
  
  path <- paste(src, fname, sep="")
  return(fromJSON(path))
  
}


load_x_standard2 <- function(location_indicators, agg_factor=8, standardize=T){
  
  caPr <- load_prism_pcs2()
  caPr.disc <- aggregate(caPr, fact=agg_factor)
  
  x_1 <- caPr.disc[[1]][]
  x_1 <- x_1[][!is.na(x_1[])]
  x_2 <- caPr.disc[[2]][]
  x_2 <- x_2[][!is.na(x_2[])]
  if (standardize){
    mu_1 <- mean(x_1[location_indicators])
    sd_1 <- sd(x_1[location_indicators])
    mu_2 <- mean(x_2[location_indicators])
    sd_2 <- sd(x_2[location_indicators])
  } else{
    mu_1 <- 0
    mu_2 <- 0
    sd_1 <- 1
    sd_2 <- 1
  }
  
  x_1_std <- (x_1 - mu_1)/sd_1
  x_2_std <- (x_2 - mu_2)/sd_2
  x_std <- array(1, c(length(x_1), 1))
  x_std <- cbind(x_std, x_1_std, x_2_std)
  return(x_std)
  
}


load_x_ca2 <- function(factor=NULL){
  
  caPr <- load_prism_pcs2()
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


make_rmse_row <- function(rmses, pattern, years, model){
  row <- list(Pattern=pattern, Model=model)
  for (i in 1:length(years)){
    row[[as.character(years[i])]] <- round(rmses[i],3)
  }
  return(row)
}


summarize_params <- function(output, true_params, species, model){
  
  bc0 <- round(mean(output$samples.beta.ca[,1]), 3)
  bc0_v <- round(var(output$samples.beta.ca[,1]), 3)
  bc1 <- round(mean(output$samples.beta.ca[,2]), 3)
  bc1_v <- round(var(output$samples.beta.ca[,2]), 3)
  bc2 <- round(mean(output$samples.beta.ca[,3]), 3)
  bc2_v <- round(var(output$samples.beta.ca[,3]), 3)
  
  bco0 <- round(mean(output$samples.beta.co[,1]), 3)
  bco0_v <- round(var(output$samples.beta.co[,1]), 3)
  bco1 <- round(mean(output$samples.beta.co[,2]), 3)
  bco1_v <- round(var(output$samples.beta.co[,2]), 3)
  bco2 <- round(mean(output$samples.beta.co[,3]), 3)
  bco2_v <- round(var(output$samples.beta.co[,3]), 3)
  alpha.ca <- round(mean(output$samples.alpha.ca), 3)
  alpha.ca_v <- round(var(output$samples.alpha.ca), 3)
  alpha.co <- round(mean(output$samples.alpha.co), 3)
  alpha.co_v <- round(var(output$samples.alpha.co), 3)
  
  rows <- list()
  rows[[1]] <- make_row_multi(bc0, true_params$beta.cases[species, 1], bc0_v, 'Beta 0 (case)', species, model)
  rows[[2]] <- make_row_multi(bc1, true_params$beta.cases[species, 2], bc1_v, 'Beta 1 (case)', species, model)
  rows[[3]] <- make_row_multi(bc2, true_params$beta.cases[species, 3], bc2_v, 'Beta 2 (case)', species, model)
  
  rows[[4]] <- make_row_multi(bco0, true_params$beta.ctrls[species, 1], bco0_v, 'Beta 0 (control)', species, model)
  rows[[5]] <- make_row_multi(bco1, true_params$beta.ctrls[species, 2], bco1_v, 'Beta 1 (control)', species, model)
  rows[[6]] <- make_row_multi(bco2, true_params$beta.ctrls[species, 3], bco2_v, 'Beta 2 (control)', species, model)
  
  rows[[7]] <- make_row_multi(alpha.ca, true_params$Alpha.cases[[1]], alpha.ca_v, 'Alpha (case)', species, model)
  rows[[8]] <- make_row_multi(alpha.co, true_params$Alpha.ctrls[[1]], alpha.co_v, 'Alpha (control)', species, model)
  
  return(rows)
  
}
