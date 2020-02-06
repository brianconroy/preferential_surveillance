

rangeMhUpdate <- function(theta, w, D, phi, proposal.sd, a, b){
  
  # proposal
  theta.new <- rlnorm(1, meanlog=log(theta), sdlog=proposal.sd)
  q.new <- dlnorm(theta.new, meanlog=log(theta), sdlog=proposal.sd, log=T)
  q.old <- dlnorm(theta, meanlog=log(theta.new), sdlog=proposal.sd, log=T)
  
  # covariance matrices
  sigma.curr <- Exponential(D, range=theta, phi=phi)
  sigma.new <- Exponential(D, range=theta.new, phi=phi)
  
  # likelihoods
  loglik.curr <- dmvnorm(w, sigma=sigma.curr, log=T)
  loglik.new <- dmvnorm(w, sigma=sigma.new, log=T)
  like.diff <- loglik.new - loglik.curr
  
  # priors
  prior.curr <- dgamma(theta, shape=a, scale=b, log=T)
  prior.new <- dgamma(theta.new, shape=a, scale=b, log=T)
  prior.diff <- prior.new - prior.curr
  
  out <- list()
  acceptance <- exp(like.diff + prior.diff + q.old - q.new)
  if(runif(1) <= acceptance) {
    out$theta <- theta.new
    out$accept <- 1
  } else { 
    out$theta <- theta
    out$accept <- 0
  }
  
  return(out)
  
}


logit <- function(x){log(x/(1-x))}


expit <- function(x){exp(x)/(1+exp(x))}


logPost <- function(y, x, beta, offset){
  
  
  prob1 <- expit(x %*% beta + offset)
  like <- sum(dbinom(y, 1, prob1, log=TRUE))
  prior <- sum(dnorm(beta, 0, 10, log=TRUE))
  return(like + prior)
  
  
}
