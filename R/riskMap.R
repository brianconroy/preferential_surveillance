
# generate posterior risk samples at high resolution
calc_posterior_risk_ds <- function(output, x, w_ds){

  n.samp <- nrow(output$samples.beta.ca)
  n.cell <- ncol(w_ds)
  risk_samp <- array(NA, c(n.samp, n.cell))

  for (i in 1:n.samp){

    beta_ca <- output$samples.beta.ca[i,]
    beta_co <- output$samples.beta.co[i,]
    alpha_ca <- output$samples.alpha.ca[i]
    alpha_co <- output$samples.alpha.co[i]
    w <- w_ds[i]

    lodds.i <- x %*% beta_ca + alpha_ca * w - x %*% beta_co - alpha_co * w
    risk_samp[i,] <- t(calc_risk(lodds.i))

  }

  return(risk_samp)

}


calc_significance <- function(samples.risk, r, threshold){
  
  pvals <- sapply(1:ncol(samples.risk), function(i){
    r_i <- samples.risk[,i]
    sum(r_i > threshold)/length(r_i)
  })
  
  # overlay p values on raster
  rp <- r
  rp[][!is.na(rp[])] <- pvals
  
  point_est <- sapply(1:ncol(samples.risk), function(i){
    r_i <- samples.risk[,i]
    mean(r_i)
  })
  
  # overlay posterior means onto raster
  r_point_est <- r
  r_point_est[][!is.na(r_point_est[])] <- point_est

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


calc_risk <- function(lodds){
  
  return(exp(lodds)/(1 + exp(lodds)))
  
}
