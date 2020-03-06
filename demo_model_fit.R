
############################
# Hamiltonian MCMC to sample
# from a bivariate normal
############################


mu <- c(0, 0)
sigma <- matrix(c(1, 0.8, 0.8, 1), nrow=2)
sigma.inv <- solve(sigma)


# step size
delta <- 0.3
nSamples <- 1000
L <- 20

# potential energy
# -log of the target density 
U <- function(x){
  return(0.5 * t(x) %*% sigma.inv %*% x)
}

# gradient of U
dU <- function(x){
  return(t(x) %*% sigma.inv)
}

# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}

# initial state 
x <- array(NA, c(nSamples, 2))
x0 <- c(0, 6)
x[1,] <- x0

for (i in 1:(nSamples-1)){
  
  # sample random momentum
  p0 <- matrix(rnorm(2))
  
  # simulate Hamiltonian dynamics
  # first 1/2 step of momentum
  xcurr <- matrix(x[i,])
  pStar <- p0 - t(0.5 * delta * dU(xcurr))
  
  # first full step for position
  xStar <- xcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - t(delta * dU(xStar))
    
    # position
    xStar <- xStar + delta*pStar
  }
 
  # last half step
  pStar <- pStar - t(0.5 * delta * dU(xStar))
  
  
  # evaluate energies
  U0 <- U(x[i,])
  UStar <- U(xStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (runif(1, 0, 1) < alpha){
    x[i+1,] <- xStar
  } else {
    x[i+1,] <- x[i,]
  }
  
}

plot(x = x[,1], y=x[,2], type='l')
points(x=x[,1], y=x[,2])
lines(x=x[1:50,1], y=x[1:50,2], type='l', col='2')
points(x=x[1:50,1], y=x[1:50,2], col='2')

colMeans(x)
var(x[,1])
var(x[,2])
cov(x[,1], x[,2])


######################
# Hamiltonian MCMC for 
# Poisson regression
######################


Beta <- c(2, 1)
N <- 100
X <- cbind(rep(1, N), rnorm(N))
rates <- exp(X %*% Beta)
Y <- sapply(rates, function(x){rpois(n=1, x)})

prior_mu <- rep(0, length(Beta))
prior_var <- rep(1000, length(Beta))
prior_sigma <- diag(prior_var)
prior_sigma_inv <- solve(prior_sigma)


# potential energy
U <- function(b){
  
  # likelihood
  logd <- 0
  lin_preds <- X %*% b
  for (i in 1:N){
    logd <- logd + dpois(Y[i], lambda=exp(lin_preds[i]), log=TRUE)
  }
  
  # prior
  for (i in 1:length(b)){
    logd <- logd + dnorm(b[i], mean=prior_mu[i], sd=sqrt(prior_var[i]), log=T)
  }
  
  return(-logd)
}


# gradient of U
dU <- function(b){
  
  grad <- array(0, c(length(Beta), 1))
  lin_preds <- X %*% b
  for (i in 1:length(Y)){
    for (j in 1:length(b)){
      grad[j] <- grad[j] + X[i,j]*Y[i] - X[i,j]*exp(lin_preds[i,])
    }
  }
  
  # prior contribution
  grad <- grad + t(t(b) %*% prior_sigma_inv)
  
  return(-grad)
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


# step size
delta <- 0.0005
nSamples <- 20000
L <- 20


# initial state 
samples.b <- array(NA, c(nSamples, length(Beta)))
b0 <- c(0, 1)
samples.b[1,] <- b0


for (i in 1:(nSamples-1)){
  
  # sample random momentum
  p0 <- matrix(rnorm(2))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(samples.b[i,])
  pStar <- p0 - 0.5 * delta * dU(bcurr)
  
  # first full step for position
  bStar <- bcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(bStar)
    
    # position
    bStar <- bStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(bStar)
  
  
  # evaluate energies
  U0 <- U(samples.b[i,])
  UStar <- U(bStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (runif(1, 0, 1) < alpha){
    samples.b[i+1,] <- bStar
  } else {
    samples.b[i+1,] <- samples.b[i,]
  }
  
}

plot(samples.b[,1], type='l'); abline(h=Beta[1], col=2)
plot(samples.b[,2], type='l'); abline(h=Beta[2], col=2)
betaHat <- colMeans(samples.b)
print(betaHat)


######################
# Hamiltonian MCMC for 
# Logistic regression
######################


logit <- function(x){
  return(log(x) - log(1-x))
}


expit <- function(x){
  return(1/(1 + exp(-x)))
}


Beta <- c(-1, 2)
N <- 100
X <- cbind(rep(1, N), rnorm(N))
probs <- expit(X %*% Beta)
Y <- sapply(probs, function(x){rbinom(n=1, size=1, x)})


prior_mu <- rep(0, length(Beta))
prior_var <- rep(1000, length(Beta))
prior_sigma <- diag(prior_var)
prior_sigma_inv <- solve(prior_sigma)


# potential energy
U <- function(b){
  
  # likelihood
  logd <- 0
  lin_preds <- X %*% b
  for (i in 1:N){
    logd <- logd + dbinom(Y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  for (i in 1:length(b)){
    logd <- logd + dnorm(b[i], mean=prior_mu[i], sd=sqrt(prior_var[i]), log=T)
  }
  
  return(-logd)
}


# gradient of U
dU <- function(b){
  
  grad <- array(0, c(length(Beta), 1))
  lin_preds <- X %*% b
  for (i in 1:length(Y)){
    for (j in 1:length(b)){
      grad[j] <- grad[j] + X[i,j]*Y[i] - X[i,j]*expit(lin_preds[i,])
    }
  }
  
  # prior contribution
  grad <- grad + t(t(b) %*% prior_sigma_inv)
  
  return(-grad)
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


# step size
delta <- 0.005
nSamples <- 20000
L <- 20


# initial state 
samples.b <- array(NA, c(nSamples, length(Beta)))
b0 <- c(0, 1)
samples.b[1,] <- b0


for (i in 1:(nSamples-1)){
  
  # sample random momentum
  p0 <- matrix(rnorm(2))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(samples.b[i,])
  pStar <- p0 - 0.5 * delta * dU(bcurr)
  
  # first full step for position
  bStar <- bcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(bStar)
    
    # position
    bStar <- bStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(bStar)
  
  
  # evaluate energies
  U0 <- U(samples.b[i,])
  UStar <- U(bStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (runif(1, 0, 1) < alpha){
    samples.b[i+1,] <- bStar
  } else {
    samples.b[i+1,] <- samples.b[i,]
  }
  
}

plot(samples.b[,1], type='l'); abline(h=Beta[1], col=2)
plot(samples.b[,2], type='l'); abline(h=Beta[2], col=2)
betaHat <- colMeans(samples.b)
print(betaHat)


############################
# Hamiltonian MCMC to sample
# from Gaussian process
# random effects
############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 5
Phi <- 3
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


probLogisticUpdate <- function(logit_prob, y, proposal_sd, x, beta, w){
  
  logit_prob.new <- logit_prob
  accept <- 0
  prior_mean <- x %*% beta + w
  
  for(i in 1:length(logit_prob)){
    
    logit_prob.prop <- rnorm(1, logit_prob[i], proposal_sd)
    prob.prop <- expit(logit_prob.prop)
    prob.curr <- expit(logit_prob[i])
    
    ll.prop <- dbinom(y[i], 1, prob.prop, log=TRUE)
    ll.curr <- dbinom(y[i], 1, prob.curr, log=TRUE)
    
    prior_mean_i <- prior_mean[i,]
    pr.prop <- dnorm(logit_prob.prop, prior_mean_i, 2, log=TRUE)
    pr.curr <- dnorm(logit_prob[i], prior_mean_i, 2, log=TRUE)
    
    acceptance <- exp(ll.prop + pr.prop - ll.curr - pr.curr)
    if(runif(1) <= acceptance) {
      logit_prob.new[i] <- logit_prob.prop
      accept <- accept + 1
    }
    
  }
  
  out <- list()
  out$logit_prob <- logit_prob.new
  out$accept <- accept
  return(out)
  
}


subdiv_matrix <- function(dmat, sigma, i, n){
  
  top <- sort(dmat[i,])[0:(n+1)]
  ids <- as.numeric(names(top))
  return(ids)
  
}


subdiv_ids <- list()
for (i in 1:length(w)){
  subdiv_ids[[i]] <- subdiv_matrix(d, sigma, i, 10)
}


subdiv_sigmas <- list()
for (i in 1:length(w)){
  ids_i <- subdiv_ids[[i]]
  subdiv_sigmas[[i]] <- sigma[ids_i, ids_i]
}


#################
## HMC functions
#################


# potential energy of a single w
# U <- function(y, mu, w, i, prior.mean, prior.var){
#   
#   mu.i <- mu[i]
#   w.i <- w[i]
#   y.i <- y[i]
#   
#   # likelihood
#   lin_pred <- mu.i + w.i
#   log_lik <- dnorm(y.i, mean = lin_pred, log=T)
#   
#   # prior
#   sigma_sub <- subdiv_sigmas[[i]]
#   ids_sub <- subdiv_ids[[i]]
#   w_sub <- w[subdiv_ids[[i]]]
#   i_char <- as.character(i)
#   i_not <- as.character(ids_sub[!ids_sub == i])
#   
#   sigma12 <- matrix(sigma_sub[i_char, i_not], nrow=1)
#   sigma22 <- sigma[i_not, i_not]
#   sigma22.inv <- solve(sigma22)
#   w_minus <- w[i_not]
#   prior_mean <- sigma12 %*% sigma22.inv %*% w_minus
#   
#   # prior variance
#   sigma11 <- sigma_sub[i_char, i_char]
#   sigma21 <- matrix(sigma_sub[i_not, i_char], ncol=1)
#   prior_var <- sigma11 - sigma12 %*% sigma22.inv %*% sigma21
#   
#   prior <- dnorm(w.i, prior_mean, sqrt(prior_var), log=T)
#   return(-(log_lik + prior))
#   
# }
# 
# nrgk <- c()
# for (i in 1:192){
#   nrgk <- c(nrgk, U(logit_prob, mu, w, i))
# }
# 
# # negate
# dU <- function(w, prior.mean, prior.var){
#   return(-(w - prior.mean)/prior.var)
# }

# 


# update all w's together

# potential energy
U <- function(y, w, mu, sigma){
  
  lin_preds <- mu + w
  return(-dmvnorm(y, lin_preds, sigma, log=T))
  
}

sigma.inv <- solve(sigma)

# gradient of potential energy
dU <- function(y, w, mu, sigma.inv){
  
  lin_preds <- mu + w
  return(t(-t(lin_preds) %*% sigma.inv))
  
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y, w, mu, sigma){
  
    
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
    
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(y, wcurr, mu, sigma.inv)
    
  # first full step for position
  wStar <- wcurr + delta*pStar
    
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y, wStar, mu, sigma.inv)
      
    # position
    wStar <- wStar + delta*pStar
  }
    
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y, wStar, mu, sigma.inv)
  
  # evaluate energies
  U0 <- U(y, wcurr, mu, sigma)
  UStar <- U(y, wStar, mu, sigma)
    
  K0 <- K(p0)
  KStar <- K(pStar)
    
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
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
  return(out)
  
}


prob <- runif(N, 0, 1)
logit_prob <- logit(prob)
x.l <- locs$x.scaled
y.l <- locs$status
mu <- x.l %*% beta.samp

delta <- 0.005
L <- 10

nrgk <- 0
for (i in 1:100){
  out <- wHmcUpdate(logit_prob, W, mu, sigma)
  nrgk <- nrgk + out$accept
}


n.sample <- 15000
N <- nrow(d)

beta.l <- beta.samp
w.i <- W
prob <- runif(N, 0, 1)
logit_prob <- logit(prob)

proposal.sd.beta.l <- 0.2
proposal.sd.w <- 0.5
proposal_sd.prob <- 0.5

samples.beta.l <- array(NA, c(n.sample, length(beta.samp)))
samples.w <- array(NA, c(n.sample, N))
samples.logit <- array(NA, c(n.sample, length(logit_prob)))

accept <- c(0, 0, 0)

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  # sample from logit prob
  logit.out <- probLogisticUpdate(logit_prob, y.l, proposal_sd.prob, x.l, beta.l, w.i)
  logit_prob <- logit.out$logit_prob
  
  ## Sample from w
  sigma.i <- sigma
  mu.i <- x.l %*% beta.l
  w.out.i <- wHmcUpdate(logit_prob, w.i, mu.i, sigma.i)
  w.i <- w.out.i$w
  w.i <- w.i - mean(w.i)
  
  accept[1] <- accept[1] + logit.out$accept
  accept[2] <- accept[2] + w.out.i$accept
  samples.logit[i,] <- logit_prob
  samples.w[i,] <- t(w.i)
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat)
abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}


############################
# Hamiltonian MCMC to sample
# from Gaussian process
# random effects. Logistic
# model update
############################

library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 6
Phi <- 4
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
# set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-0.5, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W) #, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


logit <- function(x){
  return(log(x) - log(1-x))
}


expit <- function(x){
  return(1/(1 + exp(-x)))
}


# potential energy
U <- function(y, w, mu, sigma){
  
  # likelihood
  logd <- 0
  lin_preds <- mu + w
  for (i in 1:N){
    logd <- logd + dbinom(y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


# gradient of potential energy
dU <- function(y, w, mu, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- mu + w
  
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y[i] - expit(lin_preds[i,])
  }
  
  grad <- grad + t(-t(w) %*% sigma.inv)
  return(-grad)
  
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y, w, mu, sigma, sigma.inv){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(y, wcurr, mu, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y, wStar, mu, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y, wStar, mu, sigma.inv)
  
  # evaluate energies
  U0 <- U(y, wcurr, mu, sigma)
  UStar <- U(y, wStar, mu, sigma)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
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
  return(out)
  
}

x.l <- locs$x.scaled
y.l <- locs$status
mu <- x.l %*% beta.samp

n.sample <- 8000
delta <- 0.1
L <- 12
accept <- c(0, 0, 0)

beta.l <- beta.samp
w.i <- rnorm(length(cells.all)) # W
phi.i <- Phi
theta.i <- Theta
N <- nrow(d)

prior.phi <- c(3, 6)
prior.theta <- c(2.2, 2.2)

samples.w <- array(NA, c(n.sample, N))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
proposal.sd.theta <- 0.4

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  ## Sample from w
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  sigma.inv.i <- solve(sigma.i)
  mu.i <- x.l %*% beta.l
  w.out.i <- wHmcUpdate(y.l, w.i, mu.i, sigma.i, sigma.inv.i)
  w.i <- w.out.i$w

  ## sample from Theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  theta.i <- Theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  samples.phi[i,] <- phi.i#
  phi.i <- Phi
  
  accept[2] <- accept[2] + w.out.i$accept
  accept[3] <- accept[3] + theta.out$accept
  samples.w[i,] <- t(w.i)
  samples.theta[i,] <- theta.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

accept <- accept/n.sample
print(accept)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

plot(samples.theta, type='l'); abline(h=Theta, col=2)
print(mean(samples.theta))

plot(samples.phi, type='l'); abline(h=Phi, col=2)
print(mean(samples.phi))


############################
# Hamiltonian MCMC to sample
# from Gaussian process
# random effects. 
# Poisson counts
############################

library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 5
Phi <- 3
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)


#### Simulate locations
beta.samp <- c(1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


#### Disease covariate surface
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]


simCountsGp <- function(r, loc.stats, beta, beta.samp, w, global.center=F, seed=NULL){
  
  
  # extract covariate values at survey locations
  k <- length(loc.stats$cells)
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- extract(r[[i]], loc.stats$coords)
    if (global.center){
      vals <- values(r[[i]])[!is.na(values(r[[i]]))]
      mu <- mean(vals)
      sdev <- sd(vals) 
    } else {
      mu <- mean(x.i)
      sdev <- sd(x.i)
    }
    x.i.standard <- (x.i - mu)/sdev
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }
  
  
  output <- list()

  
  # gaussian process
  w <- w[as.logical(loc.stats$status)]
  
  
  # simulate counts
  rates <- exp(x.standardised %*% beta + w)
  if (!is.null(seed)){ set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  output$y <- counts
  output$x.standardised <- x.standardised
  output$x <- x
  output$p <- ncol(x)
  output$w <- w
  return(output)
  
  
}


#### Simulate counts given locations
count.data <- simCountsGp(cov.disc, locs, beta.case, beta.samp, W, seed=42)


# potential energy
U <- function(y, w, mu, sigma){
  
  # likelihood
  logd <- 0
  lin_preds <- mu + w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


# gradient of potential energy
dU <- function(y, w, mu, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- mu + w
  
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y[i] - exp(lin_preds[i,])
  }
  
  grad <- grad + t(-t(w) %*% sigma.inv)
  return(-grad)
  
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y, w, mu, sigma, sigma.inv){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(y, wcurr, mu, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y, wStar, mu, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y, wStar, mu, sigma.inv)
  
  # evaluate energies
  U0 <- U(y, wcurr, mu, sigma)
  UStar <- U(y, wStar, mu, sigma)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
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
  return(out)
  
}


n.sample <- 8000
delta <- 0.1
L <- 20
accept <- c(0, 0, 0)

X.c <- count.data$x.standardised
Y <- count.data$y
N.y <- length(Y)
W.sub <- rnorm(sum(locs$status)) # W[as.logical(locs$status)]
sigma.sub <- Sigma[as.logical(locs$status), as.logical(locs$status)]
d.sub <- d[as.logical(locs$status), as.logical(locs$status)]

samples.w <- array(NA, c(n.sample, length(Y)))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
proposal.sd.theta <- 0.4

prior.phi <- c(3, 6)
prior.theta <- c(2.2, 2.2)

theta.i <- Theta
phi.i <- Phi
w.i <- W.sub

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for (i in 1:n.sample){
  
  ## sample from w
  sigma.i <- Exponential(d.sub, range=theta.i, phi=phi.i)
  sigma.inv.i <- solve(sigma.i)
  mu.i <- X.c %*% beta.case
  w.out.i <- wHmcUpdate(Y, w.i, mu.i, sigma.i, sigma.inv.i)
  w.i <- w.out.i$w
  
  ## sample from Theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d.sub, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta

  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N.y/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  samples.phi[i,] <- phi.i

  accept[2] <- accept[2] + w.out.i$accept
  accept[3] <- accept[3] + theta.out$accept
  samples.w[i,] <- t(w.i)
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

w.hat <- colMeans(samples.w)
w.true <- W[as.logical(locs$status)]
plot(x=w.true, y=w.hat); abline(0, 1, col=2)
plot(x=w.true, y=100*(W.sub-w.hat)/W.sub); abline(h=0, col='2')
summary(100*(w.true-w.hat)/w.true)
summary(abs(100*(w.true-w.hat)/w.true))

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=w.true[j*i], col='2')
}

plot(samples.theta, type='l'); abline(h=Theta, col='2')
print(mean(samples.theta))

plot(samples.phi, type='l'); abline(h=Phi, col='2')
print(mean(samples.phi))


############################
# Hamiltonian MCMC to sample
# from case only preferential
# sampling model
############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 5
Phi <- 5
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
# set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)


#### Simulate locations
beta.samp <- c(-2, 1)
loc.disc <- caWc.disc[[c(1)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W)#, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


#### Disease covariate surface
alpha <- 1
beta.case <- c(2, 1)
cov.disc <- caWc.disc[[c(1)]]


#### Simulate counts given locations
count.data <- simConditionalGp(cov.disc, locs, beta.case, beta.samp, alpha, seed=42)


# potential energy
U <- function(y.c, y.l, w, x.c, x.l, beta.c, beta.l, alpha, sigma, loc.stats){
  
  logd <- 0
  
  # likelihood: counts
  loc.pred <- x.l %*% beta.l + w
  offset <- alpha * expit(loc.pred)
  offset.c <- offset[as.logical(loc.stats$status)]
  rates.c <- exp(x.c %*% beta.c + offset.c)
  
  for (i in 1:length(y.c)){
    logd <- logd + dpois(y.c[i], lambda=rates.c[i], log=T)
  }
  
  # likelihood: locations
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc.pred[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


# potential energy gradient
dU <- function(y.c, y.l, w, x.c, x.l, beta.c, beta.l, sigma.inv, alpha, loc.stats){
  
  grad <- array(0, c(length(w), 1))
  lin.count <- x.c %*% beta.c
  lin.loc <- x.l %*% beta.l

  # count contribution
  ids <- c()
  for (i in 1:length(loc.stats$status)){
    if (as.logical(loc.stats$status[i])){
      ids <- c(ids, i)
    }
  }
  
  for (i in 1:length(ids)){
    id.i <- ids[i]
    f.w <- lin.count[i] + alpha * expit(lin.loc[id.i] + w[id.i])
    f.w.prime <- alpha * exp(lin.loc[id.i] + w[id.i]) / (1 + exp(lin.loc[id.i] + w[id.i]))^2
    grad[id.i] <- grad[id.i] + y.c[i] * f.w.prime - exp(f.w) * f.w.prime
  }
  
  # location contribution # WRONG!!!!
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(lin.loc[i] + w[i])
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y.c, y.l, w, x.c, x.l, beta.c, beta.l, 
                       sigma, sigma.inv, alpha, loc.stats){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(y.c, y.l, wcurr, x.c, x.l, beta.c, 
                                 beta.l, sigma.inv, alpha, loc.stats)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y.c, y.l, wStar, x.c, x.l, beta.c,
                                beta.l, sigma.inv, alpha, loc.stats)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y.c, y.l, wStar, x.c, x.l, beta.c,
                                    beta.l, sigma.inv, alpha, loc.stats)
  
  # evaluate energies
  U0 <- U(y.c, y.l, wcurr, x.c, x.l, beta.c, beta.l, alpha, sigma, loc.stats)
  UStar <- U(y.c, y.l, wStar, x.c, x.l, beta.c, beta.l, alpha, sigma, loc.stats)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
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
  return(out)
  
}


X.c <- count.data$x.standardised
X.l <- locs$x.scaled
Y.c <- count.data$y
Y.l <- locs$status
N.w <- length(locs$status)
if (ncol(X.l) > 1){
  X.l.sub <- X.l[as.logical(locs$status),]
} else{
  X.l.sub <- matrix(X.l[as.logical(locs$status),])
}

alpha.i <- alpha # runif(1, 0, 4)
theta.i <- runif(1, 3, 7)
phi.i <- runif(1, 3, 9)
w.i <- rnorm(length(cells.all))
beta.c <- beta.case # rnorm(ncol(X.c), 0, 2)
beta.l <- rnorm(ncol(X.l), 0, 2)
p.c <- length(beta.c)
p.l <- length(beta.l)

n.sample <- 3000
accept <- rep(0, 5)

samples.w <- array(NA, c(n.sample, length(Y.l)))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.alpha <- array(NA, c(n.sample, 1))
samples.beta.c <- array(NA, c(n.sample, p.c))
samples.beta.l <- array(NA, c(n.sample, p.l))

delta <- 0.1
L <- 12

proposal.sd.beta.c <- 0.04
proposal.sd.beta.l <- 0.3
proposal.sd.theta <- 0.3
proposal.sd.alpha <- 0.07

prior.phi <- c(3, 12)
prior.theta <- c(3, 3)
prior.mean.beta <- rep(0, p.c)
prior.var.beta <- rep(1000, p.c)

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for (i in 1:n.sample){
  
  ## sample from w
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  sigma.inv.i <- solve(sigma.i)
  w.out.i <- wHmcUpdate(Y.c, Y.l, w.i, X.c, X.l, beta.c, beta.l, 
                        sigma.i, sigma.inv.i, alpha.i, locs)
  w.i <- w.out.i$w

  ## sample from theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N.w/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  samples.phi[i,] <- phi.i
  
  ## sample from alpha
  # w.i.sub <- w.i[as.logical(locs$status)]
  # loc.pred <- expit(X.l.sub %*% beta.l + w.i.sub)
  # alpha.out <- alphaPoissonUpdate(X.c, Y.c, alpha.i, loc.pred, beta.c, proposal.sd.alpha)
  # alpha.i <- alpha.out$alpha
  
  ## sample from beta (disease)
  # offset.beta.c <- alpha.i * loc.pred
  # beta.out.c <- betaPoissonUpdate(X.c, Y.c, beta.c, offset.beta.c,
  #                                 prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
  # beta.c <- beta.out.c$beta
  
  ## sample from beta (location)
  # beta.out.l <- betaLocationalUpdate(X.l, Y.l, X.l.sub, X.c, Y.c, beta.l,
  #                                   beta.c, alpha.i, proposal.sd.beta.l, w.i, w.i.sub)
  beta.out.l <- betaLogisticUpdate(X.l, Y.l, beta.l, proposal.sd.beta.l, w.i)
  beta.l <- beta.out.l$beta
  
  accept[1] <- accept[1] + w.out.i$accept
  accept[2] <- accept[2] + theta.out$accept
  # accept[3] <- accept[3] + alpha.out$accept
  # accept[4] <- accept[4] + beta.out.c$accept
  accept[5] <- accept[5] + beta.out.l$accept
  
  samples.w[i,] <- t(w.i)
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  samples.alpha[i,] <- alpha.i
  samples.beta.c[i,] <- beta.c
  samples.beta.l[i,] <- beta.l
  
  ## tune the proposal standard deviations
  k <- i/100
  if(ceiling(k)==floor(k)){
    # proposal.sd.beta.l <- tuners.sd.1(accept[5], i, proposal.sd.beta.l, 50, 60)
    # proposal.sd.beta.c <- tuners.sd.1(accept[4], i, proposal.sd.beta.c, 40, 50)
    # proposal.sd.alpha <- tuners.sd.1(accept[3], i, proposal.sd.alpha, 50, 60)
  }

  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}


accept <- accept/n.sample
print(accept)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
plot(x=W, y=100*(W-w.hat)/W)

summary(100*(W-w.hat)/W)
summary(abs(100*(W-w.hat)/W))
summary(100*(W[abs(W)>1]-w.hat[abs(W)>1])/W[abs(W)>1])

par(mfrow=c(4, 4))
j <- 2
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

plot(samples.theta, type='l'); abline(h=Theta, col='2')
hist(samples.theta); abline(v=Theta, col='2')
print(mean(samples.theta))

plot(samples.phi, type='l'); abline(h=Phi, col='2')
hist(samples.phi); abline(v=Phi, col='2')
print(mean(samples.phi))

plot(samples.alpha, type='l'); abline(h=alpha, col='2')
hist(samples.alpha); abline(v=alpha, col='2')

par(mfrow=c(1,2))
plot(samples.beta.c[,1], type='l'); abline(h=beta.case[1], col='2')
plot(samples.beta.c[,2], type='l'); abline(h=beta.case[2], col='2')

plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col='2')
plot(samples.beta.l[,2], type='l'); abline(h=beta.samp[2], col='2')



############################
# Hamiltonian MCMC to sample
# from Gaussian process
# random effects. Logistic
# model update
############################

library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 5
Phi <- 5
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 2)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W) #, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


logit <- function(x){
  return(log(x) - log(1-x))
}


expit <- function(x){
  return(1/(1 + exp(-x)))
}


# potential energy
U <- function(y, w, mu, sigma){
  
  # likelihood
  logd <- 0
  lin_preds <- mu + w
  for (i in 1:N){
    logd <- logd + dbinom(y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


# gradient of potential energy
dU <- function(y, w, mu, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- mu + w
  
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y[i] - expit(lin_preds[i,])
  }
  
  grad <- grad + t(-t(w) %*% sigma.inv)
  return(-grad)
  
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y, w, mu, sigma, sigma.inv){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(y, wcurr, mu, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y, wStar, mu, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y, wStar, mu, sigma.inv)
  
  # evaluate energies
  U0 <- U(y, wcurr, mu, sigma)
  UStar <- U(y, wStar, mu, sigma)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
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
  return(out)
  
}

x.l <- locs$x.scaled
y.l <- locs$status

n.sample <- 3000
delta <- 0.4
L <- 15
accept <- c(0, 0, 0)

p.l <- length(beta.samp)
beta.l <- rnorm(p.l, 0, 2)
w.i <- rnorm(length(cells.all))
phi.i <- runif(1, 2, 8)
theta.i <- runif(1, 2, 8)
N <- nrow(d)

prior.phi <- c(3, 6)
prior.theta <- c(2.2, 2.2)
prior.mean.beta <- rep(0, p.l)
prior.var.beta <- rep(1000, p.l)

samples.w <- array(NA, c(n.sample, N))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.beta.l <- array(NA, c(n.sample, p.l))

proposal.sd.theta <- 0.4
proposal.sd.beta.l <- 0.4

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  ## Sample from w
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  sigma.inv.i <- solve(sigma.i)
  mu.i <- x.l %*% beta.l
  w.out.i <- wHmcUpdate(y.l, w.i, mu.i, sigma.i, sigma.inv.i)
  w.i <- w.out.i$w
  
  ## sample from theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  theta.i <- Theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  phi.i <- Phi
  
  ## sample from beta
  beta.out.l <- betaLogisticUpdate(x.l, y.l, beta.l, proposal.sd.beta.l, w.i)
  beta.l <- beta.out.l$beta
  beta.l <- beta.samp
  
  accept[1] <- accept[1] + beta.out.l$accept
  accept[2] <- accept[2] + w.out.i$accept
  accept[3] <- accept[3] + theta.out$accept
  samples.w[i,] <- t(w.i)
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  samples.beta.l[i,] <- beta.l
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

accept <- accept/n.sample
print(accept)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

plot(samples.theta, type='l'); abline(h=Theta, col=2)
print(mean(samples.theta))

plot(samples.phi, type='l'); abline(h=Phi, col=2)
print(mean(samples.phi))

plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col=2)
hist(samples.beta.l[,1]); abline(v=beta.samp[1], col=2)
print(mean(samples.beta.l[,1]))

plot(samples.beta.l[,2], type='l'); abline(h=beta.samp[2], col=2)
hist(samples.beta.l[,2]); abline(v=beta.samp[2], col=2)
print(mean(samples.beta.l[,2]))


############################
# Hamiltonian MCMC to sample
# from Gaussian process
# random effects. Joint 
# beta location / w update
############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 5
Phi <- 5
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 2)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


logit <- function(x){
  return(log(x) - log(1-x))
}


expit <- function(x){
  return(1/(1 + exp(-x)))
}


# potential energy
U <- function(y, x, w_b, n_w, sigma){
  
  w <- w_b[1:n_w]
  beta <- w_b[(n_w + 1):nrow(w_b)]
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + w
  for (i in 1:N){
    logd <- logd + dbinom(y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  logd <- logd + dmvnorm(beta, rep(0, length(beta)), diag(rep(5, length(beta))), log=T)
  
  return(-logd)
  
}


# gradient of potential energy
dU <- function(y, x, w_b, n_w, sigma.inv){
  
  w <- w_b[1:n_w]
  beta <- w_b[(n_w + 1):nrow(w_b)]
  
  grad <- array(0, c(length(w) + length(beta), 1))
  lin_preds <- x %*% beta + w
  n_w <- length(w)
  
  for (i in 1:n_w){
    grad[i] <- grad[i] + y[i] - expit(lin_preds[i,])
  }
  
  for (j in 1:length(beta)){
    for (i in 1:n_w){
      grad[n_w + j] <- grad[n_w + j] + x[i, j] * (y[i] -  expit(lin_preds[i,]))
    }
  }
  
  prior_beta_inv <- diag(rep(1/5, length(beta)))
  
  grad[1:n_w] <- grad[1:n_w] + t(-t(w) %*% sigma.inv)
  grad[(n_w+1):nrow(grad)] <- grad[(n_w+1):nrow(grad)] + t(-t(beta) %*% prior_beta_inv)
  
  return(-grad)
  
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y, w, x, beta, sigma, sigma.inv){
  
  
  # sample random momentum
  n_w <- length(w)
  n_b <- length(beta)
  p0 <- matrix(rnorm(n_w + n_b))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  betacurr <- matrix(beta)
  wbcurr <- rbind(wcurr, betacurr)
  pStar <- p0 - 0.5 * delta * dU(y, x, wbcurr, n_w, sigma.inv)
  
  # first full step for position
  wbStar <- wbcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y, x, wbStar, n_w, sigma.inv)
    
    # position
    wbStar <- wbStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y, x, wbStar, n_w, sigma.inv)
  
  # evaluate energies
  U0 <- U(y, x, wbcurr, n_w, sigma)
  UStar <- U(y, x, wbStar, n_w, sigma)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (runif(1, 0, 1) < alpha){
    wnext <- wbStar[1:n_w]
    bnext <- wbStar[(n_w+1):nrow(wbStar)]
    accept <- 1
  } else {
    wnext <- wcurr
    bnext <- betacurr
    accept <- 0
  }
  
  out <- list()
  out$w <- wnext
  out$beta <- bnext
  out$accept <- accept
  return(out)
  
}

x.l <- locs$x.scaled
y.l <- locs$status

n.sample <- 3000
delta <- 0.3
L <- 5
accept <- c(0, 0, 0)

p.l <- length(beta.samp)
beta.l <- rnorm(p.l, 0, 2)
w.i <- rnorm(length(cells.all))
phi.i <- runif(1, 2, 8)
theta.i <- runif(1, 2, 8)
N <- nrow(d)
n_w <- length(w.i)

prior.phi <- c(3, 6)
prior.theta <- c(2.2, 2.2)
prior.mean.beta <- rep(0, p.l)
prior.var.beta <- rep(1000, p.l)

samples.w <- array(NA, c(n.sample, N))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.beta.l <- array(NA, c(n.sample, p.l))

proposal.sd.theta <- 0.4
proposal.sd.beta.l <- 0.4

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  ## Sample from w and beta
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  sigma.inv.i <- solve(sigma.i)
  wb.out.i <- wHmcUpdate(y.l, w.i, x.l, beta.l, sigma.i, sigma.inv.i)
  w.i <- wb.out.i$w
  beta.l <- wb.out.i$beta
  
  ## sample from theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  
  ## sample from beta
  # beta.out.l <- betaLogisticUpdate(x.l, y.l, beta.l, proposal.sd.beta.l, w.i)
  # beta.l <- beta.out.l$beta
  
  # accept[1] <- accept[1] + beta.out.l$accept
  accept[2] <- accept[2] + wb.out.i$accept
  accept[3] <- accept[3] + theta.out$accept
  samples.w[i,] <- t(w.i)
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  samples.beta.l[i,] <- beta.l
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

accept <- accept/n.sample
print(accept)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

plot(samples.theta, type='l'); abline(h=Theta, col=2)
print(mean(samples.theta))

plot(samples.phi, type='l'); abline(h=Phi, col=2)
print(mean(samples.phi))

plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col=2)
hist(samples.beta.l[,1]); abline(v=beta.samp[1], col=2)
print(mean(samples.beta.l[,1]))

plot(samples.beta.l[,2], type='l'); abline(h=beta.samp[2], col=2)
hist(samples.beta.l[,2]); abline(v=beta.samp[2], col=2)
print(mean(samples.beta.l[,2]))



############################
# Hamiltonian MCMC to sample
# from Gaussian process
# random effects. 
# Separate hamiltonian updates
############################

library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 8
Phi <- 5
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-0.5, 1)
loc.disc <- caWc.disc[[c(1)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


logit <- function(x){
  return(log(x) - log(1-x))
}


expit <- function(x){
  return(1/(1 + exp(-x)))
}


# potential energy
U <- function(y, w, mu, sigma){
  
  # likelihood
  logd <- 0
  lin_preds <- mu + w
  for (i in 1:N){
    logd <- logd + dbinom(y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


Ubeta <- function(y, w, x, beta){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + w
  for (i in 1:N){
    logd <- logd + dbinom(y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  prior_beta_mean <- rep(0, length(beta))
  prior_beta_var <- diag(rep(5, length(beta)))
  logd <- logd + dmvnorm(as.numeric(beta), prior_beta_mean, prior_beta_var, log=T)
  
  return(-logd)
  
}




# gradient of potential energy
dU <- function(y, w, mu, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- mu + w
  
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y[i] - expit(lin_preds[i,])
  }
  
  grad <- grad + t(-t(w) %*% sigma.inv)
  return(-grad)
  
}


# gradient of potential energy
dU.beta <- function(y, w, x, beta){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- x %*% beta + w
  
  for (j in 1:length(beta)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x[i, j] * (y[i] - expit(lin_preds[i,]))
    }
  }
  
  prior_var_inv <- diag(rep(1/5, length(beta)))
  grad <- grad + t(-t(beta) %*% prior_var_inv)
  return(-grad)
  
}


# kinetic energy
K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y, w, mu, sigma, sigma.inv){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(y, wcurr, mu, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y, wStar, mu, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y, wStar, mu, sigma.inv)
  
  # evaluate energies
  U0 <- U(y, wcurr, mu, sigma)
  UStar <- U(y, wStar, mu, sigma)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
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
  return(out)
  
}


betaHmcUpdate <- function(y, w, x, beta){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta)
  pStar <- p0 - 0.5 * delta_b * dU.beta(y, w, x, bcurr)
  
  # first full step for position
  bStar <- bcurr + delta_b*pStar
  
  # full steps
  for (jL in 1:c(L_b-1)){
    # momentum
    pStar <- pStar - delta_b * dU.beta(y, w, x, bStar)
    
    # position
    bStar <- bStar + delta_b*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_b * dU.beta(y, w, x, bStar)
  
  # evaluate energies
  U0 <- Ubeta(y, w, x, bcurr)
  UStar <- Ubeta(y, w, x, bStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
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
  return(out)
  
}


x.l <- locs$x.scaled
y.l <- locs$status

n.sample <- 3000
delta <- 0.3
L <- 12

delta_b <- 0.3
L_b <- 12
accept <- c(0, 0, 0)

beta.l <- rnorm(2)
w.i <- rnorm(length(cells.all))
phi.i <- Phi
theta.i <- Theta
N <- nrow(d)

prior.phi <- c(3, 6)
prior.theta <- c(2.2, 2.2)

samples.w <- array(NA, c(n.sample, N))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.beta.l <- array(NA, c(n.sample, length(beta.l)))
proposal.sd.theta <- 0.4

progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  ## sample from w
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  sigma.inv.i <- solve(sigma.i)
  mu.i <- x.l %*% beta.l
  w.out.i <- wHmcUpdate(y.l, w.i, mu.i, sigma.i, sigma.inv.i)
  w.i <- w.out.i$w
  w.i <- w.i - mean(w.i)
  # 
  # ## sample from theta
  # theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  # theta.i <- theta.out$theta
  # 
  # ## sample from phi
  # R.i <- sigma.i/phi.i
  # phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  # 
  ## sample from beta
  beta.out.l <- betaHmcUpdate(y.l, w.i, x.l, beta.l)
  beta.l <- beta.out.l$beta
  
  accept[1] <- accept[1] + beta.out.l$accept
  accept[2] <- accept[2] + w.out.i$accept
  #accept[3] <- accept[3] + theta.out$accept
  samples.w[i,] <- t(w.i)
  #samples.phi[i,] <- phi.i
  #samples.theta[i,] <- theta.i
  samples.beta.l[i,] <- beta.l
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

accept <- accept/n.sample
print(accept)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

plot(samples.theta, type='l'); abline(h=Theta, col=2)
print(mean(samples.theta))

plot(samples.phi, type='l'); abline(h=Phi, col=2)
print(mean(samples.phi))

plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col=2)
plot(samples.beta.l[,2], type='l'); abline(h=beta.samp[2], col=2)
