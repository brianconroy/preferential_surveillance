

calc_posterior_risk <- function(output, samples_int){
  
  samples.beta.ca <- output$samples.beta.ca
  samples.beta.co <- output$samples.beta.co
  samples.alpha.ca <- output$samples.alpha.ca
  samples.alpha.co <- output$samples.alpha.co
  X_high <- load_x_ca()
  n.samples <- nrow(samples_int)
  samples.risk <- c()
  
  samples.l <- lapply(1:n.samples, function(i){
    w_i <- samples_int[i,]
    beta.ca_i <- samples.beta.ca[i,]
    beta.co_i <- samples.beta.co[i,]
    alpha.ca_i <- samples.alpha.ca[i]
    alpha.co_i <- samples.alpha.co[i]
    lodds_i <- X_high %*% beta.ca_i + alpha.ca_i * w_i - X_high %*% beta.co_i - alpha.co_i * w_i
    t(calc_risk(lodds_i))
  })
  samples.risk <- do.call(rbind, samples.l)
  
  return(samples.risk)
  
}


calc_significance <- function(samples.risk, r, threshold){
  
  pvals <- sapply(1:ncol(samples.risk), function(i){
    r_i <- samples.risk[,i]
    sum(r_i > threshold)/length(r_i)
  })
  
  rp <- overlay(pvals, r)
  
  point_est <- sapply(1:ncol(samples.risk), function(i){
    r_i <- samples.risk[,i]
    mean(r_i)
  })
  
  r_point_est <- overlay(point_est, r)
  
  inds_95 <- 1 * (pvals > 0.95)
  inds_50 <- 1 * (pvals > 0.5)
  inds_25 <- 1 * (pvals > 0.25)
  # normalize color scales
  if (sum(inds_95) == 0){inds_95[1]=3}
  if (sum(inds_50) == 0){inds_95[2]=2}
  if (sum(inds_25) == 0){inds_95[3]=1}
  
  inds <- inds_95 + inds_50 + inds_25
  r_inds <- overlay(inds, r)
  r_inds_95 <- overlay(inds_95, r)
  
  return(list(
    r_point_est=r_point_est,
    r_fracs=rp,
    r_inds=r_inds,
    r_inds_95=r_inds_95
  ))
  
}
