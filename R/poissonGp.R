

poissonGp <- function(x, y, d, n.sample, burnin, L_w=20, L_b=22, proposal.sd.theta=0.3,
                      w_initial=NULL, beta_initial=NULL, theta_initial=NULL, phi_initial=NULL,
                      prior_phi=c(4, 12), prior_theta=c(2.5, 2.5)){
  
  
  # initial values
  N.w <- ncol(d)
  p.b <- ncol(x)
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(beta_initial)){
    beta.i <- rnorm(p.b)
  } else {
    beta.i <- beta_initial
  } 
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(phi_initial)){
    phi.i <- runif(1, 3, 5)
  } else {
    phi.i <- phi_initial
  }
  
  # storage
  accept <- rep(0, 3)
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, N.w))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.beta <- array(NA, c(n.keep, p.b))
  
  # dual averaging quantities
  w_tuning <- initialize_tuning(m=700, target=0.75)
  b_tuning <- initialize_tuning(m=2000, target=0.75)
  
  deltas_w <- c()
  deltas_b <- c()
  
  # priors
  prior.mean.beta <- rep(0, p.b)
  prior.var.beta <- rep(1000, p.b)
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateBench(x, y, beta.i, w.i, sigma.i, sigma.inv.i, w_tuning$delta_curr, L_w)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N.w/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])
    
    ## sample from beta
    beta.out.i <- caseHmcUpdateBench(y, w.i, x, beta.i, b_tuning$delta_curr, L_b)
    beta.i <- beta.out.i$beta
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.beta[j,] <- beta.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept[1] <- accept[1] + w.out.i$accept
      accept[2] <- accept[2] + theta.out$accept
      accept[3] <- accept[3] + beta.out.i$accept
      
    }
    
    w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
    deltas_w <- c(deltas_w, w_tuning$delta_curr)
    
    b_tuning <- update_tuning(b_tuning, beta.out.i$a, i, beta.out.i$accept)
    deltas_b <- c(deltas_b, b_tuning$delta_curr)
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  
  output <- list()
  output$accept <- accept
  output$samples.beta <- samples.beta
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$deltas_b <- deltas_b
  output$n.sample <- n.sample
  output$burnin <- burnin
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  
  return(output)
  
  
}


burnin_poisson_gp <- function(output_poisson, n.burn){
  
  n.curr <- nrow(output_poisson$samples.w)
  i.start <- n.burn + 1
  output_poisson$n.sample <- output_poisson$n.sample - n.burn
  output_poisson$burnin <- output_poisson$burnin + n.burn
  output_poisson$samples.w <- output_poisson$samples.w[i.start:n.curr,]
  output_poisson$samples.phi <- output_poisson$samples.phi[i.start:n.curr]
  output_poisson$samples.beta <- output_poisson$samples.beta[i.start:n.curr,]
  output_poisson$samples.theta <- output_poisson$samples.theta[i.start:n.curr]
  return(output_poisson)
  
}


continue_poisson_gp <- function(data, x, y, output_poisson, n.sample, d.sub){
  
  # get initial values
  n.sample.old <- nrow(output_poisson$samples.w)
  w_initial <- output_poisson$samples.w[n.sample.old,]
  theta_initial <- output_poisson$samples.theta[n.sample.old]
  phi_initial <- output_poisson$samples.phi[n.sample.old]
  beta_initial <- output_poisson$samples.beta[n.sample.old,]
  
  # get tuning parameters
  delta_w <- tail(output_poisson$deltas_w, 1)
  L_w <- 8
  deltas_beta <- tail(output_poisson$deltas_b, 1)
  L_b <- 8
  proposal.sd.theta <- output_poisson$proposal.sd.theta
  prior_phi <- output_poisson$prior_phi
  prior_theta <- output_poisson$prior_theta
  
  # run fit
  more_output_poisson <- poissonGp(x, y, d.sub,
                                   n.sample=n.sample, burnin=0, proposal.sd.theta=proposal.sd.theta,
                                   L_w=L_w, L_b=L_b,
                                   beta_initial=beta_initial, w_initial=w_initial, 
                                   phi_initial=phi_initial, theta_initial=theta_initial,
                                   prior_phi=prior_phi, prior_theta=prior_theta)
  
  # merge samples
  new_output_poisson <- output_poisson
  new_output_poisson$samples.phi <- c(new_output_poisson$samples.phi, more_output_poisson$samples.phi)
  new_output_poisson$samples.theta <- c(new_output_poisson$samples.theta, more_output_poisson$samples.theta)
  new_output_poisson$samples.w <- rbind(new_output_poisson$samples.w, more_output_poisson$samples.w)
  new_output_poisson$samples.beta <- rbind(new_output_poisson$samples.beta, more_output_poisson$samples.beta)
  new_output_poisson$deltas_w <- c(new_output_poisson$deltas_w, more_output_poisson$deltas_w)
  new_output_poisson$deltas_b <- c(new_output_poisson$deltas_b, more_output_poisson$deltas_b)
  new_output_poisson$n.sample <- new_output_poisson$n.sample + n.sample
  new_output_poisson$accept <- (more_output_poisson$accept * n.sample + output_poisson$accept * output_poisson$n.sample)/(n.sample + output_poisson$n.sample)
  
  return(new_output_poisson)
  
}


Uw_bench <- function(x, y, beta, w, sigma){
  
  logd <- 0
  
  # likelihood: counts
  count_pred <- x %*% beta + w
  rates <- exp(count_pred)
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


Ubeta_bench <- function(y, w, x, beta){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + w
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


dU_w_bench <- function(x, y, beta, w, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  
  # count contribution
  lin.count <- x %*% beta
  for (i in 1:length(y)){
    grad[i] <- grad[i] + (y[i] - exp(lin.count[i] + w[i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


dU_beta_bench <- function(y, w, x, beta){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- x %*% beta + w
  
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


wHmcUpdateBench <- function(x, y, beta, w, sigma, sigma.inv, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_bench(x, y, beta, wcurr, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_bench(x, y, beta, wStar, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_bench(x, y, beta, wStar, sigma.inv)
  
  # evaluate energies
  U0 <- Uw_bench(x, y, beta, wcurr, sigma)
  UStar <- Uw_bench(x, y, beta, wStar, sigma)
  
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


caseHmcUpdateBench <- function(y, w, x, beta, delta_c, L_c){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta)
  pStar <- p0 - 0.5 * delta_c * dU_beta_bench(y, w, x, bcurr)
  
  # first full step for position
  bStar <- bcurr + delta_c*pStar
  
  # full steps
  for (jL in 1:c(L_c-1)){
    # momentum
    pStar <- pStar - delta_c * dU_beta_bench(y, w, x, bStar)
    
    # position
    bStar <- bStar + delta_c*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_c * dU_beta_bench(y, w, x, bStar)
  
  # evaluate energies
  U0 <- Ubeta_bench(y, w, x, bcurr)
  UStar <- Ubeta_bench(y, w, x, bStar)
  
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


krigeW <- function(output, d, ids){
  
  n.new <- nrow(d) - length(ids)
  nsamp.old <- nrow(output$samples.w)
  samples.new <- array(NA, c(nrow(output$samples.w), n.new))
  cat("generating", nsamp.old, "kriged samples\n")
  progressBar <- txtProgressBar(style=3)
  percentage.points<-round((1:100/100)*nsamp.old)
  
  for (i in 1:nsamp.old){
    
    theta.i <- output$samples.theta[i]
    phi.i <- output$samples.phi[i]
    covmat.i <- Exponential(d, range=theta.i, phi=phi.i)
    
    omega.11 <- covmat.i[-ids, -ids]
    omega.12 <- covmat.i[-ids, ids]
    omega.21 <- covmat.i[ids, -ids]
    omega.22 <- covmat.i[ids, ids]
    
    w.i <- output$samples.w[i,]
    e.cond <- omega.12 %*% solve(omega.22) %*% w.i
    var.cond <- omega.11 - omega.12 %*% solve(omega.22) %*% omega.21
    w.pred <- mvrnorm(n=1, mu=e.cond, var.cond)
    samples.new[i,] <- w.pred
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/nsamp.old)
    }
    
  }
  
  out <- list(
    mu.new=colMeans(samples.new),
    samples.new=samples.new
  )
  
  return(out)
  
}
