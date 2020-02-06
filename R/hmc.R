

Uw <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma, loc.stats){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- x.c %*% beta.c + alpha * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.c)){
    logd <- logd + dpois(y.c[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


Ubeta <- function(y, w, x, beta, alpha, offset){
  
  # likelihood
  logd <- 0
  lin_preds <- offset + x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  prior_beta_mean <- rep(0, length(beta))
  if (length(beta) > 1){
    prior_beta_var <- diag(rep(100, length(beta)))
  } else {
    prior_beta_var <- matrix(100)
  }
  
  logd <- logd + dmvnorm(as.numeric(beta), prior_beta_mean, prior_beta_var, log=T)
  
  return(-logd)
  
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


dU_w <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma.inv, loc.stats){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(w[i])
  }
  
  # count contribution
  lin.count <- x.c %*% beta.c
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha * (y.c[i] - exp(lin.count[i] + alpha * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


dU_beta <- function(y, w, x, beta, alpha, offset){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- offset + x %*% beta + alpha * w
  
  for (j in 1:length(beta)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x[i, j] * (y[i] - exp(lin_preds[i,]))
    }
  }
  
  if (length(beta) > 1){
    prior_var_inv <- diag(rep(1/100, length(beta)))
  } else{
    prior_var_inv <- 1/100
  }
  
  grad <- grad + t(-t(beta) %*% prior_var_inv)
  return(-grad)
  
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


K <- function(p){
  return(t(p) %*% p/2)
}


wHmcUpdate <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma, sigma.inv, loc.stats, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w(y.l, x.c, y.c, alpha, beta.c, wcurr, sigma.inv, loc.stats)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w(y.l, x.c, y.c, alpha, beta.c, wStar, sigma.inv, loc.stats)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w(y.l, x.c, y.c, alpha, beta.c, wStar, sigma.inv, loc.stats)
  
  # evaluate energies
  U0 <- Uw(y.l, x.c, y.c, alpha, beta.c, wcurr, sigma, loc.stats)
  UStar <- Uw(y.l, x.c, y.c, alpha, beta.c, wStar, sigma, loc.stats)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    alpha <- 0
  }
  
  if (runif(1, 0, 1) < alpha){
    wnext <- wStar
    accept <- 1
  } else {
    wnext <- wcurr
    accept <- 0
  }
  
  out <- list()
  out$w <- wnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}


betaHmcUpdate <- function(y, w, x, beta, alpha, delta_c, L_c, offset){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta)
  pStar <- p0 - 0.5 * delta_c * dU_beta(y, w, x, bcurr, alpha, offset)
  
  # first full step for position
  bStar <- bcurr + delta_c*pStar
  
  # full steps
  for (jL in 1:c(L_c-1)){
    # momentum
    pStar <- pStar - delta_c * dU_beta(y, w, x, bStar, alpha, offset)
    
    # position
    bStar <- bStar + delta_c*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_c * dU_beta(y, w, x, bStar, alpha, offset)
  
  # evaluate energies
  U0 <- Ubeta(y, w, x, bcurr, alpha, offset)
  UStar <- Ubeta(y, w, x, bStar, alpha, offset)
  
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
  out$beta <- bnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}


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


Uw_cc <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma, loc.stats, offset){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: case counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- offset + x.c %*% beta.ca + alpha.ca * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=rates[i], log=T)
  }
  
  # likelihood: control counts
  count_pred <- offset + x.c %*% beta.co + alpha.co * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.co)){
    logd <- logd + dpois(y.co[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_cc <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma.inv, loc.stats, offset){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(w[i])
  }
  
  # case contribution
  lin.count <- offset + x.c %*% beta.ca
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.ca * (y.ca[i] - exp(lin.count[i] + alpha.ca * w[id.i]))
  }
  
  # control contribution
  lin.count <- offset + x.c %*% beta.co
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.co * (y.co[i] - exp(lin.count[i] + alpha.co * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateCC <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma, sigma.inv, loc.stats, delta, L, offset){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wcurr, sigma.inv, loc.stats, offset)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats, offset)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats, offset)
  
  # evaluate energies
  U0 <- Uw_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wcurr, sigma, loc.stats, offset)
  UStar <- Uw_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma, loc.stats, offset)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    alpha <- 0
  }
  
  if (runif(1, 0, 1) < alpha){
    wnext <- wStar
    accept <- 1
  } else {
    wnext <- wcurr
    accept <- 0
  }
  
  out <- list()
  out$w <- wnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}
