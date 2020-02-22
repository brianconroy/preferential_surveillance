

betaLocHmcUpdate <- function(y.l, w, x.l, beta.loc, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta.loc)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta.loc)
  pStar <- p0 - 0.5 * delta * dU_beta_loc(y.l, w, x.l, beta.loc)
  
  # first full step for position
  bStar <- bcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_beta_loc(y.l, w, x.l, bStar)
    
    # position
    bStar <- bStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_beta_loc(y.l, w, x.l, bStar)
  
  # evaluate energies
  U0 <- U_beta_loc(y.l, w, x.l, bcurr)
  UStar <- U_beta_loc(y.l, w, x.l, bStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    alpha <- 0
  }
  
  if (runif(1, 0, 1) < alpha){
    bnext <- bStar
    accept <- 1
  } else {
    bnext <- bcurr
    accept <- 0
  }
  
  out <- list()
  out$b <- bnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}


dU_beta_loc <- function(y.l, w, x.l, beta.loc){
  
  grad <- array(0, c(length(beta.loc), 1))
  
  # location contribution
  loc_pred <- x.l %*% beta.loc + w
  for (j in 1:length(beta.loc)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x.l[i,j] * (y.l[i] - expit(loc_pred[i]))
    }
  }
  
  # prior contribution
  prior_var_inv <- diag(rep(1/100, length(beta.loc)))
  grad <- grad + t(-t(beta.loc) %*% prior_var_inv)
  return(-grad)
  
}


U_beta_loc <- function(y.l, w, x.l, beta.loc){
  
  # likelihood
  logd <- 0
  lin_preds <- x.l %*% beta.loc + w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  prior_beta_mean <- rep(0, length(beta.loc))
  prior_beta_var <- diag(rep(100, length(beta.loc)))
  
  logd <- logd + dmvnorm(as.numeric(beta.loc), prior_beta_mean, prior_beta_var, log=T)
  
  return(-logd)
  
}
