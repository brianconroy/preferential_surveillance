

logisticGp <- function(y, x, d, n.sample, burnin, L_beta, L_w, proposal.sd.theta=0.3,
                       w_initial=NULL, theta_initial=NULL, phi_initial=NULL, beta_initial=NULL,
                       prior_phi, prior_theta){
  
  
  # initial values
  N.w <- ncol(d)
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(phi_initial)){
    phi.i <- runif(1, 10, 15)
  } else {
    phi.i <- phi_initial
  }
  if (is.null(beta_initial)){
    beta.i <- rep(0, ncol(x))
  } else {
    beta.i <- beta_initial
  }
  
  # storage
  accept <- list(w=0, theta=0, beta=0)
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, N.w))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.beta <- array(NA, c(n.keep, ncol(x)))

  # dual averaging quantities
  w_tuning <- initialize_tuning(m=700, target=0.75)
  deltas_w <- c()
  beta_tuning <- initialize_tuning(m=700, target=0.75)
  deltas_beta <- c()
  
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from beta
    beta.out.i <- betaLocHmcUpdate(y, w.i, x, beta.i, beta_tuning$delta_curr, L_beta)
    beta.i <- beta.out.i$b
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateLogitCov(y, w.i, x, beta.i, sigma.i, sigma.inv.i, w_tuning$delta_curr, L_w)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N.w/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.theta[j,] <- theta.i
      samples.beta[j,] <- beta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept$theta <- accept$theta + theta.out$accept
      accept$beta <- accept$beta + beta.out.i$accept
      accept$w <- accept$w + w.out.i$accept
      
    }
    
    w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
    deltas_w <- c(deltas_w, w_tuning$delta_curr)
    
    beta_tuning <- update_tuning(beta_tuning, beta.out.i$a, i, beta.out.i$accept)
    deltas_beta <- c(deltas_beta, beta_tuning$delta_curr)
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  for (n in names(accept)){
    accept[[n]] <- accept[[n]]/n.keep
  }

  output <- list()
  output$accept <- accept
  output$samples.theta <- samples.theta
  output$samples.beta <- samples.beta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$deltas_beta  <- deltas_beta
  output$burnin <- burnin
  output$n.sample <- n.sample
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$L_w <- L_w
  output$L_beta <- L_beta
  output$proposal.sd.theta <- proposal.sd.theta

  return(output)
  
  
}


burnin_logistic_gp_cov <- function(output, n.burn){
  
  n.curr <- output$n.sample - output$burnin
  i.start <- n.burn + 1
  output$burnin <- output$burnin + n.burn
  output$samples.w <- output$samples.w[i.start:n.curr,]
  output$samples.phi <- output$samples.phi[i.start:n.curr]
  output$samples.beta <- output$samples.beta[i.start:n.curr,]
  output$samples.theta <- output$samples.theta[i.start:n.curr]
  return(output)
  
}


continue_logistic_gp_cov <- function(data, output, n.sample){
  
  # get initial values
  n.sample.old <- nrow(output$samples.w)
  w_initial <- output$samples.w[n.sample.old,]
  theta_initial <- output$samples.theta[n.sample.old]
  phi_initial <- output$samples.phi[n.sample.old]
  beta_initial <- output$samples.beta[n.sample.old,]
  
  # get tuning parameters
  delta_w <- tail(output$deltas_w, 1)
  L_w <- output$L_w
  deltas_beta <- tail(output$deltas_beta, 1)
  L_beta <- output$L_beta
  proposal.sd.theta <- output$proposal.sd.theta
  prior_phi <- output$prior_phi
  prior_theta <- output$prior_theta
  
  more_output <- logisticGpCov(data$loc$status, data$loc$x.scaled, 
                               d, n.sample, burnin=0, L_beta=L_beta, L_w=L_w, proposal.sd.theta=proposal.sd.theta,
                               w_initial=w_initial, theta_initial=theta_initial, phi_initial=phi_initial,
                               beta_initial=beta_initial,
                               prior_phi=prior_phi, prior_theta=prior_theta)
  
  # merge samples
  new_output <- output
  new_output$samples.phi <- c(new_output$samples.phi, more_output$samples.phi)
  new_output$samples.theta <- c(new_output$samples.theta, more_output$samples.theta)
  new_output$samples.w <- rbind(new_output$samples.w, more_output$samples.w)
  new_output$samples.beta <- rbind(new_output$samples.beta, more_output$samples.beta)
  new_output$deltas_w <- c(new_output$deltas_w, more_output$deltas_w)
  new_output$deltas_beta <- c(new_output$deltas_beta, more_output$deltas_beta)
  new_output$n.sample <- new_output$n.sample + n.sample
  for (n in names(new_output$accept)){
    new_output$accept[[n]] <- (new_output$n.sample*new_output$accept[[n]] + n.sample*more_output$accept[[n]])/(new_output$n.sample + n.sample)
  }
  
  return(new_output)
  
}


Uw_logit_cov <- function(y, w, x, beta, sigma){
  
  logd <- 0
  
  # likelihood
  probs <- expit(x%*%beta + w)
  for (i in 1:length(y)){
    logd <- logd + dbinom(y[i], size=1, prob=probs[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_logit_cov <- function(y, w, x, beta, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  
  # likelihood contribution
  lin_pred <- x %*% beta + w
  for (i in 1:length(y)){
    grad[i] <- grad[i] + (y[i] - expit(lin_pred[i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateLogitCov <- function(y, w, x, beta, sigma, sigma.inv, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_logit_cov(y, wcurr, x, beta, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_logit_cov(y, wStar, x, beta, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_logit_cov(y, wStar, x, beta, sigma.inv)
  
  # evaluate energies
  U0 <- Uw_logit_cov(y, wcurr, x, beta, sigma)
  UStar <- Uw_logit_cov(y, wStar, x, beta, sigma)
  
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
