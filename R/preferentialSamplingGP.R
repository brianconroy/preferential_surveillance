

burnin_after <- function(output, n.burn){
  
  n.curr <- output$n.sample - output$burnin
  i.start <- n.burn + 1
  output$burnin <- output$burnin + n.burn
  if (length(dim(output$samples.alpha.ca)) == 2 | class(output$samples.alpha.ca) == 'numeric'){
    output$samples.alpha.ca <- output$samples.alpha.ca[i.start:n.curr]
    output$samples.alpha.co <- output$samples.alpha.co[i.start:n.curr]
    output$samples.beta.ca <- output$samples.beta.ca[i.start:n.curr,]
    output$samples.beta.co <- output$samples.beta.co[i.start:n.curr,]
  } else {
    n_s <- dim(output$samples.alpha.ca)[1]
    new.aca <- array(NA, c(n_s,n.curr-i.start+1,1))
    new.aco <- array(NA, c(n_s,n.curr-i.start+1,1))
    for (j in 1:n_s){
      new.aca[j,,] <- matrix(output$samples.alpha.ca[j,i.start:n.curr,])
      new.aco[j,,] <- matrix(output$samples.alpha.co[j,i.start:n.curr,])
    }
    output$samples.alpha.ca <- new.aca
    output$samples.alpha.co <- new.aco
    output$samples.beta.ca <- output$samples.beta.ca[,i.start:n.curr,]
    output$samples.beta.co <- output$samples.beta.co[,i.start:n.curr,]
  }
  output$samples.w <- output$samples.w[i.start:n.curr,]
  output$samples.phi <- output$samples.phi[i.start:n.curr]
  output$samples.theta <- output$samples.theta[i.start:n.curr]
  return(output)
  
}


#' continueMCMC
#' 
#' continues running an MCMC chain from the output of prefSampleGpCC
#'
#' @param output (list) output of prefSampleGpCC
#' @param n.sample (numeric) number of additional samples to generate
#'
#' @return
#' @export
#'
#' @examples
continueMCMC <- function(data, D, output, n.sample){
  
  # get initial values
  n.sample.old <- nrow(output$samples.beta.ca)
  beta_ca_initial <- output$samples.beta.ca[n.sample.old,]
  beta_co_initial <- output$samples.beta.co[n.sample.old,]
  alpha_ca_initial <- output$samples.alpha.ca[n.sample.old]
  alpha_co_initial <- output$samples.alpha.co[n.sample.old]

  w_initial <- output$samples.w[n.sample.old,]
  theta_initial <- output$samples.theta[n.sample.old]
  phi_initial <- output$samples.phi[n.sample.old]
  
  # get tuning parameters
  delta_w <- tail(output$deltas_w, 1)
  delta_aca <- tail(output$deltas_aca, 1)
  delta_aco <- tail(output$deltas_aco, 1)
  delta_ca <- tail(output$deltas_ca, 1)
  delta_co <- tail(output$deltas_co, 1)
  L_w <- output$L_w
  L_ca <- output$L_ca
  L_co <- output$L_co
  L_a_ca <- output$L_a_ca 
  L_a_co <- output$L_a_co
  proposal.sd.theta <- output$proposal.sd.theta
  
  # get priors
  prior_phi <- output$prior_phi
  prior_theta <- output$prior_theta
  prior_alpha_ca_var <- output$prior_alpha_ca_var
  prior_alpha_co_var <- output$prior_alpha_co_var
  more_output <- prefSampleGpCC(data, D, n.sample, burnin=0, 
                                  L_w=L_w, L_ca=L_ca, L_co=L_co, L_a_ca=L_a_ca, L_a_co=L_a_co,
                                  proposal.sd.theta=proposal.sd.theta,
                                  self_tune_w=FALSE, self_tune_aca=FALSE, self_tune_aco=FALSE, self_tune_ca=FALSE, self_tune_co=FALSE,
                                  delta_w=delta_w, delta_aca=delta_aca, delta_aco=delta_aco, delta_ca=delta_ca, delta_co=delta_co, 
                                  beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, 
                                  alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                                  theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                                  prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, 
                                  prior_alpha_co_var=prior_alpha_ca_var)
  
  # combine outputs
  new_output <- output
  new_output$samples.alpha.ca <- c(new_output$samples.alpha.ca, more_output$samples.alpha.ca)
  new_output$samples.alpha.co <- c(new_output$samples.alpha.co, more_output$samples.alpha.co)
  new_output$samples.beta.ca <- rbind(new_output$samples.beta.ca, more_output$samples.beta.ca)
  new_output$samples.beta.co <- rbind(new_output$samples.beta.co, more_output$samples.beta.co)
  new_output$samples.phi <- c(new_output$samples.phi, more_output$samples.phi)
  new_output$samples.theta <- c(new_output$samples.theta, more_output$samples.theta)
  new_output$samples.w <- rbind(new_output$samples.w, more_output$samples.w)
  new_output$deltas_aca <- c(new_output$deltas_aca, more_output$deltas_aca)
  new_output$deltas_aco <- c(new_output$deltas_aco, more_output$deltas_aco)
  new_output$deltas_co <- c(new_output$deltas_co, more_output$deltas_co)
  new_output$deltas_ca <- c(new_output$deltas_ca, more_output$deltas_ca)
  new_output$deltas_w <- c(new_output$deltas_w, more_output$deltas_w)
  new_output$n.sample <- new_output$n.sample + n.sample
  new_output$accept <- (new_output$n.sample * new_output$accept + n.sample * more_output$accept)/(new_output$n.sample + n.sample)
  
  return(new_output)
  
}


#' prefSampleGp
#' 
#' Fits the preferential sampling model with locations
#' driven by gaussian field (no intercept) via
#' Hamiltonian Monte Carlo samplers. Tunes step sizes
#' with dual averaging (Hoffman and Gelman, 2014).
#'
#' @param data list with elements: locs, conditional
#' @param n.sample number of mcmc samples
#' @param burnin mcmc burnin
#' @param L_w hmc simulation length for w
#' @param L_c hmc simulation length beta (case)
#' @param L_a hmc simulation length for alpha
#' @param proposal.sd.theta MH-RW proposal stdev for theta
#' @param m_a self tuning iterations for alpha
#' @param m_c self tuning iterations for beta (case)
#' @param m_w self tuning iterations for w
#' @param target_a alpha target acceptance rate
#' @param target_c beta (case) target acceptance rate
#' @param target_w w target acceptance rate
#' @param self_tune_w whether to self tune gaussian process
#' @param self_tune_a whether to self tune alpha
#' @param self_tune_c whether to self tune beta
#' @param delta_w w step size
#' @param delta_a alpha step size
#' @param delta_c beta (case) step size
#' @param beta_initial optional beta starting value
#' @param alpha_initial optional alpha starting value
#'
#' @return list of mcmc samples and related quantities
#' @export
#'
#' @examples
prefSampleGp <- function(data, n.sample, burnin, L_w, L_c, L_a, proposal.sd.theta=0.3,
                         m_a=700, m_c=700, m_w=700, target_a=0.75,
                         target_c=0.75, target_w=0.75, self_tune_w=TRUE,
                         self_tune_a=TRUE, self_tune_c=TRUE,
                         delta_w=NULL, delta_a=NULL, delta_c=NULL, beta_initial=NULL, alpha_initial=NULL){
  
  
  ## setup
  count.data <- data$conditional
  locs <- data$loc
  X.c <- count.data$x.standardised
  Y.c <- count.data$y
  Y.l <- locs$status
  N.w <- length(locs$status)
  
  ## starting values
  w.i <- rnorm(N.w)
  if (is.null(beta_initial)){
    mod <- glm(Y.c ~ X.c + w.i[locs$ids] - 1, family='poisson')
    beta.c <- as.numeric(coef(mod)[1:2])
  } else {
    beta.c <- beta_initial
  }
  if (is.null(alpha_initial)){
    alpha.i <- runif(1, 1, 3)
  } else {
    alpha.i <- alpha_initial
  }
  theta.i <- runif(1, 1, 4)
  phi.i <- runif(1, 1, 4)
  p.c <- length(beta.c)
  
  # storage
  accept <- rep(0, 4)
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, length(Y.l)))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.alpha <- array(NA, c(n.keep, 1))
  samples.beta.c <- array(NA, c(n.keep, length(beta.c)))
  
  if (self_tune_w){
    w_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_a){
    a_tuning <- initialize_tuning(m=m_a, target=target_a)
  } else {
    a_tuning <- list(delta_curr=delta_a)
  }
  if (self_tune_c){
    c_tuning <- initialize_tuning(m=m_c, target=target_c)
  } else {
    c_tuning <- list(delta_curr=delta_c)
  }
  
  deltas_w <- c()
  deltas_a <- c()
  deltas_c <- c()
  
  prior.phi <- c(4, 12)
  prior.theta <- c(2.5, 2.5)
  prior.mean.beta <- rep(0, p.c)
  prior.var.beta <- rep(1000, p.c)
  
  cat("Generating", n.sample, " samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdate(Y.l, X.c, Y.c, alpha.i, beta.c, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L_w)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N.w/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
    
    ## sample from beta.c
    w.i.sub <- w.i[locs$ids]
    beta.out.c <- betaHmcUpdate(Y.c, w.i[locs$ids], X.c, beta.c, alpha.i, c_tuning$delta_curr, L_c)
    beta.c <- beta.out.c$beta
    
    ## sample from alpha
    alpha.out <- alphaHmcUpdate(Y.c, w.i.sub, X.c, beta.c, alpha.i, a_tuning$delta_curr, L_a)
    alpha.i <- alpha.out$alpha
    
    if (i > burnin){
      
      j <- i - burnin
      samples.beta.c[j,] <- beta.c
      samples.alpha[j,] <- alpha.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept[1] <- accept[1] + w.out.i$accept
      accept[2] <- accept[2] + theta.out$accept
      accept[3] <- accept[3] + beta.out.c$accept
      accept[4] <- accept[4] + alpha.out$accept
      
    }
    
    if (self_tune_w) {
      w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
      deltas_w <- c(deltas_w, w_tuning$delta_curr)
    }
    if (self_tune_a){
      a_tuning <- update_tuning(a_tuning, alpha.out$a, i, alpha.out$accept)
      deltas_a <- c(deltas_a, a_tuning$delta_curr)
    }
    if (self_tune_c){
      c_tuning <- update_tuning(c_tuning, beta.out.c$a, i, beta.out.c$accept)
      deltas_c <- c(deltas_c, c_tuning$delta_curr)
    }

    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  
  output <- list()
  output$accept <- accept
  output$samples.beta.c <- samples.beta.c
  output$samples.alpha <- samples.alpha
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$deltas_a <- deltas_a
  output$deltas_c <- deltas_c
  
  return(output)
  
}


#' prefSampleGpCC
#' 
#' uses normal priors for alpha case and alpha control
#'
#' @param data 
#' @param d 
#' @param n.sample 
#' @param burnin 
#' @param L_w 
#' @param L_ca 
#' @param L_co 
#' @param L_a_ca 
#' @param L_a_co 
#' @param proposal.sd.theta 
#' @param m_aca 
#' @param m_aco 
#' @param m_ca 
#' @param m_co 
#' @param m_w 
#' @param target_aca 
#' @param target_aco 
#' @param target_ca 
#' @param target_co 
#' @param target_w 
#' @param self_tune_w 
#' @param self_tune_aca 
#' @param self_tune_aco 
#' @param self_tune_ca 
#' @param self_tune_co 
#' @param delta_w 
#' @param delta_aca 
#' @param delta_aco 
#' @param delta_ca 
#' @param delta_co 
#' @param beta_ca_initial 
#' @param beta_co_initial 
#' @param alpha_ca_initial 
#' @param alpha_co_initial 
#' @param theta_initial 
#' @param phi_initial 
#' @param w_initial 
#' @param prior_phi 
#' @param prior_theta 
#' @param prior_alpha_ca_var 
#' @param prior_alpha_co_var 
#'
#' @return
#' @export
#'
#' @examples
prefSampleGpCC <- function(data, d, n.sample, burnin, 
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=0.3,
                           m_aca=2000, m_aco=2000, m_ca=700, m_co=700, m_w=700, 
                           target_aca=0.75, target_aco=0.75, target_ca=0.75, target_co=0.75, target_w=0.75, 
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                           delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                           beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                           theta_initial=NULL, phi_initial=NULL, w_initial=NULL,
                           prior_phi, prior_theta, prior_alpha_ca_var, prior_alpha_co_var,
                           offset=1){
  
  
  ## setup
  case.data <- data$case.data
  ctrl.data <- data$ctrl.data
  locs <- data$loc
  X.c <- case.data$x.standardised
  Y.ca <- case.data$y
  Y.co <- ctrl.data$y
  Y.l <- locs$status
  N.w <- length(locs$status)
  log_offset <- log(offset)
  
  
  ## starting values
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(beta_ca_initial)){
    beta.ca <- rnorm(ncol(X.c))
  } else {
    beta.ca <- beta_ca_initial
  }
  if (is.null(beta_co_initial)){
    beta.co <- rnorm(ncol(X.c))
  } else {
    beta.co <- beta_co_initial
  }
  if (is.null(alpha_ca_initial)){
    alpha.ca.i <- runif(1, 1, 3)
  } else {
    alpha.ca.i <- alpha_ca_initial
  }
  if (is.null(alpha_co_initial)){
    alpha.co.i <- runif(1, -3, -1)
  } else {
    alpha.co.i <- alpha_co_initial
  }
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(phi_initial)){
    phi.i <- runif(1, 3.5, 4.5)
  } else {
    phi.i <- phi_initial
  }
  p.c <- length(beta.ca)
  
  
  # storage
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, length(Y.l)))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.alpha.ca <- array(NA, c(n.keep, 1))
  samples.alpha.co <- array(NA, c(n.keep, 1))
  samples.beta.ca <- array(NA, c(n.keep, length(beta.ca)))
  samples.beta.co <- array(NA, c(n.keep, length(beta.co)))
  
  
  if (self_tune_w){
    w_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_aca){
    a_ca_tuning <- initialize_tuning(m=m_aca, target=target_aca)
  } else {
    a_ca_tuning <- list(delta_curr=delta_aca)
  }
  if (self_tune_aco){
    a_co_tuning <- initialize_tuning(m=m_aco, target=target_aco)
  } else {
    a_co_tuning <- list(delta_curr=delta_aco)
  }
  if (self_tune_ca){
    ca_tuning <- initialize_tuning(m=m_ca, target=target_ca)
  } else {
    ca_tuning <- list(delta_curr=delta_ca)
  }
  if (self_tune_co){
    co_tuning <- initialize_tuning(m=m_co, target=target_co)
  } else {
    co_tuning <- list(delta_curr=delta_co)
  }
  
  deltas_w <- c()
  deltas_ca <- c()
  deltas_co <- c()
  deltas_aca <- c()
  deltas_aco <- c()
  
  accept <- rep(0, 6)
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateCC(Y.l, X.c, Y.ca, alpha.ca.i, beta.ca, Y.co,
                            alpha.co.i, beta.co, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L_w, offset=log_offset)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N.w/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])
    
    ## sample from beta.case
    w.i.sub <- w.i[locs$ids]
    beta.out.ca <- betaHmcUpdate(Y.ca, w.i[locs$ids], X.c, beta.ca, alpha.ca.i, ca_tuning$delta_curr, L_ca, offset=log_offset)
    beta.ca <- beta.out.ca$beta
    
    ## sample from alpha case
    alpha.out.ca <- alphaHmcUpdate(Y.ca, w.i.sub, X.c, beta.ca, alpha.ca.i, 
                                   a_ca_tuning$delta_curr, prior_alpha_ca_mean, prior_alpha_ca_var, L_a_ca, offset=log_offset)
    alpha.ca.i <- alpha.out.ca$alpha
    
    ## sample from beta.ctrl
    beta.out.co <- betaHmcUpdate(Y.co, w.i[locs$ids], X.c, beta.co, alpha.co.i, co_tuning$delta_curr, L_co, offset=log_offset)
    beta.co <- beta.out.co$beta
    
    ## sample from alpha control
    alpha.out.co <- alphaHmcUpdate(Y.co, w.i.sub, X.c, beta.co, alpha.co.i, 
                                   a_co_tuning$delta_curr, prior_alpha_co_mean, prior_alpha_co_var, L_a_co, offset=log_offset)
    alpha.co.i <- alpha.out.co$alpha
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.beta.ca[j,] <- beta.ca
      samples.beta.co[j,] <- beta.co
      samples.alpha.ca[j,] <- alpha.ca.i
      samples.alpha.co[j,] <- alpha.co.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept[1] <- accept[1] + w.out.i$accept
      accept[2] <- accept[2] + theta.out$accept
      accept[3] <- accept[3] + beta.out.ca$accept
      accept[4] <- accept[4] + beta.out.co$accept
      accept[5] <- accept[5] + alpha.out.ca$accept
      accept[6] <- accept[6] + alpha.out.co$accept
      
    }
    
    if (self_tune_w){
      w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
      deltas_w <- c(deltas_w, w_tuning$delta_curr)
    }
    if (self_tune_aca){
      a_ca_tuning <- update_tuning(a_ca_tuning, alpha.out.ca$a, i, alpha.out.ca$accept)
      deltas_aca <- c(deltas_aca, a_ca_tuning$delta_curr)
    }
    if (self_tune_aco){
      a_co_tuning <- update_tuning(a_co_tuning, alpha.out.co$a, i, alpha.out.co$accept)
      deltas_aco <- c(deltas_aco, a_co_tuning$delta_curr)
    }
    if (self_tune_ca){
      ca_tuning <- update_tuning(ca_tuning, beta.out.ca$a, i, beta.out.ca$accept)
      deltas_ca <- c(deltas_ca, ca_tuning$delta_curr)
    }
    if (self_tune_co){
      co_tuning <- update_tuning(co_tuning, beta.out.co$a, i, beta.out.co$accept)
      deltas_co <- c(deltas_co, co_tuning$delta_curr)
    }
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  
  output <- list()
  output$accept <- accept
  output$samples.beta.ca <- samples.beta.ca
  output$samples.beta.co <- samples.beta.co
  output$samples.alpha.ca <- samples.alpha.ca
  output$samples.alpha.co <- samples.alpha.co
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$deltas_aca <- deltas_aca
  output$deltas_aco <- deltas_aco
  output$deltas_ca <- deltas_ca
  output$deltas_co <- deltas_co
  output$L_w <- L_w
  output$L_ca <- L_ca
  output$L_co <- L_co
  output$L_a_ca <- L_a_ca 
  output$L_a_co <- L_a_co
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}


prefSampleGpCC_constrained <- function(data, d, n.sample, burnin, 
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=0.3,
                           m_aca=2000, m_aco=2000, m_ca=700, m_co=700, m_w=700, 
                           target_aca=0.75, target_aco=0.75, target_ca=0.75, target_co=0.75, target_w=0.75, 
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                           delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                           beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                           theta_initial=NULL, phi_initial=NULL, w_initial=NULL,
                           prior_phi, prior_theta, prior_alpha_ca_var, prior_alpha_co_var){
  
  
  ## setup
  case.data <- data$case.data
  ctrl.data <- data$ctrl.data
  locs <- data$loc
  X.c <- case.data$x.standardised
  Y.ca <- case.data$y
  Y.co <- ctrl.data$y
  Y.l <- locs$status
  N.w <- length(locs$status)

  
  ## starting values
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(beta_ca_initial)){
    beta.ca <- rnorm(ncol(X.c))
  } else {
    beta.ca <- beta_ca_initial
  }
  if (is.null(beta_co_initial)){
    beta.co <- rnorm(ncol(X.c))
  } else {
    beta.co <- beta_co_initial
  }
  if (is.null(alpha_ca_initial)){
    alpha.ca.i <- runif(1, 1, 3)
  } else {
    alpha.ca.i <- alpha_ca_initial
  }
  if (is.null(alpha_co_initial)){
    alpha.co.i <- runif(1, -3, -1)
  } else {
    alpha.co.i <- alpha_co_initial
  }
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(phi_initial)){
    phi.i <- runif(1, 3.5, 4.5)
  } else {
    phi.i <- phi_initial
  }
  p.c <- length(beta.ca)
  
  
  # storage
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, length(Y.l)))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.alpha.ca <- array(NA, c(n.keep, 1))
  samples.alpha.co <- array(NA, c(n.keep, 1))
  samples.beta.ca <- array(NA, c(n.keep, length(beta.ca)))
  samples.beta.co <- array(NA, c(n.keep, length(beta.co)))
  
  
  if (self_tune_w){
    w_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_aca){
    a_ca_tuning <- initialize_tuning(m=m_aca, target=target_aca)
  } else {
    a_ca_tuning <- list(delta_curr=delta_aca)
  }
  if (self_tune_aco){
    a_co_tuning <- initialize_tuning(m=m_aco, target=target_aco)
  } else {
    a_co_tuning <- list(delta_curr=delta_aco)
  }
  if (self_tune_ca){
    ca_tuning <- initialize_tuning(m=m_ca, target=target_ca)
  } else {
    ca_tuning <- list(delta_curr=delta_ca)
  }
  if (self_tune_co){
    co_tuning <- initialize_tuning(m=m_co, target=target_co)
  } else {
    co_tuning <- list(delta_curr=delta_co)
  }
  
  deltas_w <- c()
  deltas_ca <- c()
  deltas_co <- c()
  deltas_aca <- c()
  deltas_aco <- c()
  
  accept <- rep(0, 6)
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateCC(Y.l, X.c, Y.ca, alpha.ca.i, beta.ca, Y.co,
                            alpha.co.i, beta.co, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L_w, offset=0)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N.w/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])
    
    ## sample from beta.case
    w.i.sub <- w.i[locs$ids]
    beta.out.ca <- betaHmcUpdate(Y.ca, w.i[locs$ids], X.c, beta.ca, alpha.ca.i, ca_tuning$delta_curr, L_ca, offset=0)
    beta.ca <- beta.out.ca$beta
    
    ## sample from beta.ctrl
    beta.out.co <- betaHmcUpdate(Y.co, w.i[locs$ids], X.c, beta.co, alpha.co.i, co_tuning$delta_curr, L_co, offset=0)
    beta.co <- beta.out.co$beta
    
    ## sample from alphas
    which_alpha <- runif(1, 0, 1)
    if (which_alpha > 0.5){
      
      # update alpha case first
      alpha.out.ca <- alphaHmcUpdate(Y.ca, w.i.sub, X.c, beta.ca, alpha.ca.i, 
                                     a_ca_tuning$delta_curr, prior_alpha_ca_mean, prior_alpha_ca_var, L_a_ca, offset=0)
      alpha.ca.i <- alpha.out.ca$alpha
      
      alpha.out.co <- alphaHmcUpdate_constrained(Y.co, w.i.sub, X.c, beta.co, alpha.co.i, constraint=alpha.ca.i, direction='less',
                                     a_co_tuning$delta_curr, prior_alpha_co_mean, prior_alpha_co_var, L_a_co, offset=0)
      alpha.co.i <- alpha.out.co$alpha

    } else {
      # update alpha ctrl first
      alpha.out.co <- alphaHmcUpdate(Y.co, w.i.sub, X.c, beta.co, alpha.co.i, 
                                     a_co_tuning$delta_curr, prior_alpha_co_mean, prior_alpha_co_var, L_a_co, offset=0)
      alpha.co.i <- alpha.out.co$alpha
      
      alpha.out.ca <- alphaHmcUpdate_constrained(Y.ca, w.i.sub, X.c, beta.ca, alpha.ca.i, constraint=alpha.co.i, direction='greater',
                                     a_ca_tuning$delta_curr, prior_alpha_ca_mean, prior_alpha_ca_var, L_a_ca, offset=0)
      alpha.ca.i <- alpha.out.ca$alpha
      
    }
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.beta.ca[j,] <- beta.ca
      samples.beta.co[j,] <- beta.co
      samples.alpha.ca[j,] <- alpha.ca.i
      samples.alpha.co[j,] <- alpha.co.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept[1] <- accept[1] + w.out.i$accept
      accept[2] <- accept[2] + theta.out$accept
      accept[3] <- accept[3] + beta.out.ca$accept
      accept[4] <- accept[4] + beta.out.co$accept
      accept[5] <- accept[5] + alpha.out.ca$accept
      accept[6] <- accept[6] + alpha.out.co$accept
      
    }
    
    if (self_tune_w){
      w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
      deltas_w <- c(deltas_w, w_tuning$delta_curr)
    }
    if (self_tune_aca){
      a_ca_tuning <- update_tuning(a_ca_tuning, alpha.out.ca$a, i, alpha.out.ca$accept)
      deltas_aca <- c(deltas_aca, a_ca_tuning$delta_curr)
    }
    if (self_tune_aco){
      a_co_tuning <- update_tuning(a_co_tuning, alpha.out.co$a, i, alpha.out.co$accept)
      deltas_aco <- c(deltas_aco, a_co_tuning$delta_curr)
    }
    if (self_tune_ca){
      ca_tuning <- update_tuning(ca_tuning, beta.out.ca$a, i, beta.out.ca$accept)
      deltas_ca <- c(deltas_ca, ca_tuning$delta_curr)
    }
    if (self_tune_co){
      co_tuning <- update_tuning(co_tuning, beta.out.co$a, i, beta.out.co$accept)
      deltas_co <- c(deltas_co, co_tuning$delta_curr)
    }
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  
  output <- list()
  output$accept <- accept
  output$samples.beta.ca <- samples.beta.ca
  output$samples.beta.co <- samples.beta.co
  output$samples.alpha.ca <- samples.alpha.ca
  output$samples.alpha.co <- samples.alpha.co
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$deltas_aca <- deltas_aca
  output$deltas_aco <- deltas_aco
  output$deltas_ca <- deltas_ca
  output$deltas_co <- deltas_co
  output$L_w <- L_w
  output$L_ca <- L_ca
  output$L_co <- L_co
  output$L_a_ca <- L_a_ca 
  output$L_a_co <- L_a_co
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}


prefSampleGpCC_truncnorm <- function(data, n.sample, burnin, 
                                  L_w, L_ca, L_co, L_a_ca, L_a_co,
                                  proposal.sd.theta=0.3,
                                  m_aca=2000, m_aco=2000, m_ca=700, m_co=700, m_w=700, 
                                  target_aca=0.75, target_aco=0.75, target_ca=0.75, target_co=0.75, target_w=0.75, 
                                  self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                                  delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                                  beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                                  theta_initial=NULL, phi_initial=NULL, w_initial=NULL,
                                  prior_phi, prior_theta, prior_alpha_ca, prior_alpha_co){
  
  
  ## setup
  case.data <- data$case.data
  ctrl.data <- data$ctrl.data
  locs <- data$loc
  X.c <- case.data$x.standardised
  Y.ca <- case.data$y
  Y.co <- ctrl.data$y
  Y.l <- locs$status
  N.w <- length(locs$status)
  
  
  ## starting values
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(beta_ca_initial)){
    beta.ca <- rnorm(ncol(X.c))
  } else {
    beta.ca <- beta_ca_initial
  }
  if (is.null(beta_co_initial)){
    beta.co <- rnorm(ncol(X.c))
  } else {
    beta.co <- beta_co_initial
  }
  if (is.null(alpha_ca_initial)){
    alpha.ca.i <- runif(1, 1, 3)
  } else {
    alpha.ca.i <- alpha_ca_initial
  }
  if (is.null(alpha_co_initial)){
    alpha.co.i <- runif(1, -3, -1)
  } else {
    alpha.co.i <- alpha_co_initial
  }
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(phi_initial)){
    phi.i <- runif(1, 3.5, 4.5)
  } else {
    phi.i <- phi_initial
  }
  p.c <- length(beta.ca)
  
  
  # storage
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, length(Y.l)))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.alpha.ca <- array(NA, c(n.keep, 1))
  samples.alpha.co <- array(NA, c(n.keep, 1))
  samples.beta.ca <- array(NA, c(n.keep, length(beta.ca)))
  samples.beta.co <- array(NA, c(n.keep, length(beta.co)))
  
  
  if (self_tune_w){
    w_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_aca){
    a_ca_tuning <- initialize_tuning(m=m_aca, target=target_aca)
  } else {
    a_ca_tuning <- list(delta_curr=delta_aca)
  }
  if (self_tune_aco){
    a_co_tuning <- initialize_tuning(m=m_aco, target=target_aco)
  } else {
    a_co_tuning <- list(delta_curr=delta_aco)
  }
  if (self_tune_ca){
    ca_tuning <- initialize_tuning(m=m_ca, target=target_ca)
  } else {
    ca_tuning <- list(delta_curr=delta_ca)
  }
  if (self_tune_co){
    co_tuning <- initialize_tuning(m=m_co, target=target_co)
  } else {
    co_tuning <- list(delta_curr=delta_co)
  }
  
  deltas_w <- c()
  deltas_ca <- c()
  deltas_co <- c()
  deltas_aca <- c()
  deltas_aco <- c()
  
  prior.mean.beta <- rep(0, p.c)
  prior.var.beta <- rep(1000, p.c)
  prior_alpha_ca_mean <- prior_alpha_ca[1]
  prior_alpha_ca_var <- prior_alpha_ca[2]
  prior_alpha_co_mean <- prior_alpha_co[1]
  prior_alpha_co_var <- prior_alpha_co[2]
  
  accept <- rep(0, 6)
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateCC(Y.l, X.c, Y.ca, alpha.ca.i, beta.ca, Y.co,
                            alpha.co.i, beta.co, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L_w)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])
    
    ## sample from beta.case
    w.i.sub <- w.i[locs$ids]
    beta.out.ca <- betaHmcUpdate(Y.ca, w.i[locs$ids], X.c, beta.ca, alpha.ca.i, ca_tuning$delta_curr, L_ca)
    beta.ca <- beta.out.ca$beta
    
    ## sample from alpha case
    alpha.out.ca <- alphaHmcUpdate_truncnorm(Y.ca, w.i.sub, X.c, beta.ca, alpha.ca.i, prior_alpha_ca_mean,
                                         prior_alpha_ca_var, a_ca_tuning$delta_curr, 
                                         bound=0, bound_type='lower', L_a_ca)
    alpha.ca.i <- alpha.out.ca$alpha
    
    ## sample from beta.ctrl
    beta.out.co <- betaHmcUpdate(Y.co, w.i[locs$ids], X.c, beta.co, alpha.co.i, co_tuning$delta_curr, L_co)
    beta.co <- beta.out.co$beta
    
    ## sample from alpha control
    alpha.out.co <- alphaHmcUpdate_truncnorm(Y.co, w.i.sub, X.c, beta.co, alpha.co.i, 
                                             prior_alpha_co_mean, prior_alpha_co_var,
                                             a_co_tuning$delta_curr, bound=0, bound_type='upper', L_a_co)
    alpha.co.i <- alpha.out.co$alpha
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.beta.ca[j,] <- beta.ca
      samples.beta.co[j,] <- beta.co
      samples.alpha.ca[j,] <- alpha.ca.i
      samples.alpha.co[j,] <- alpha.co.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept[1] <- accept[1] + w.out.i$accept
      accept[2] <- accept[2] + theta.out$accept
      accept[3] <- accept[3] + beta.out.ca$accept
      accept[4] <- accept[4] + beta.out.co$accept
      accept[5] <- accept[5] + alpha.out.ca$accept
      accept[6] <- accept[6] + alpha.out.co$accept
      
    }
    
    if (self_tune_w){
      w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
      deltas_w <- c(deltas_w, w_tuning$delta_curr)
    }
    if (self_tune_aca){
      a_ca_tuning <- update_tuning(a_ca_tuning, alpha.out.ca$a, i, alpha.out.ca$accept)
      deltas_aca <- c(deltas_aca, a_ca_tuning$delta_curr)
    }
    if (self_tune_aco){
      a_co_tuning <- update_tuning(a_co_tuning, alpha.out.co$a, i, alpha.out.co$accept)
      deltas_aco <- c(deltas_aco, a_co_tuning$delta_curr)
    }
    if (self_tune_ca){
      ca_tuning <- update_tuning(ca_tuning, beta.out.ca$a, i, beta.out.ca$accept)
      deltas_ca <- c(deltas_ca, ca_tuning$delta_curr)
    }
    if (self_tune_co){
      co_tuning <- update_tuning(co_tuning, beta.out.co$a, i, beta.out.co$accept)
      deltas_co <- c(deltas_co, co_tuning$delta_curr)
    }
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  
  output <- list()
  output$accept <- accept
  output$samples.beta.ca <- samples.beta.ca
  output$samples.beta.co <- samples.beta.co
  output$samples.alpha.ca <- samples.alpha.ca
  output$samples.alpha.co <- samples.alpha.co
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$deltas_aca <- deltas_aca
  output$deltas_aco <- deltas_aco
  output$deltas_ca <- deltas_ca
  output$deltas_co <- deltas_co
  output$L_w <- L_w
  output$L_ca <- L_ca
  output$L_co <- L_co
  output$L_a_ca <- L_a_ca 
  output$L_a_co <- L_a_co
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_mean <- prior_alpha_ca_mean
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_mean <- prior_alpha_co_mean
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}
