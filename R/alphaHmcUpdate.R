

alphaHmcUpdate <- function(y, w, x, beta, alpha, delta_a, prior_mean, prior_var, L_a, offset){
  
  
  # sample random momentum
  p0 <- rnorm(1)
  
  # simulate Hamiltonian dynamics
  acurr <- alpha
  pStar <- p0 - 0.5 * delta_a * dU_alpha(y, w, x, beta, acurr, prior_mean, prior_var, offset)
  
  # first full step for position
  aStar <- acurr + delta_a*pStar
  
  # full steps
  for (jL in 1:c(L_a-1)){
    # momentum
    pStar <- pStar - delta_a * dU_alpha(y, w, x, beta, aStar, prior_mean, prior_var, offset)
    
    # position
    aStar <- aStar + delta_a*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_a * dU_alpha(y, w, x, beta, aStar, prior_mean, prior_var, offset)
  
  # evaluate energies
  U0 <- Ualpha(y, w, x, beta, acurr, prior_mean, prior_var, offset)
  UStar <- Ualpha(y, w, x, beta, aStar, prior_mean, prior_var, offset)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  a <- min(1, exp((U0 + K0) - (UStar + KStar)))
  
  if (is.na(a)){
    a <- 0
  }
  
  if (runif(1, 0, 1) < a){
    anext <- aStar
    accept <- 1
  } else {
    anext <- acurr
    accept <- 0
  }
  
  out <- list()
  out$alpha <- anext
  out$accept <- accept
  out$a <- a
  return(out)
  
}


Ualpha <- function(y, w, x, beta, alpha, prior_mean, prior_var, offset){
  
  # likelihood
  logd <- 0
  lin_preds <- offset + x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dnorm(alpha, prior_mean, prior_var, log=T)
  
  return(-logd)
  
}


dU_alpha <- function(y, w, x, beta, alpha, prior_mean, prior_var, offset){
  
  grad <- 0
  lin_preds <- offset + x %*% beta + alpha * w
  
  for (i in 1:length(w)){
    grad <- grad + w[i] * (y[i] - exp(lin_preds[i,]))
  }
  
  prior_var_inv <- 1/prior_var
  
  grad <- grad - (alpha - prior_mean) * prior_var_inv
  return(-grad)
  
}
