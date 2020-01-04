library(fields)
library(MASS)


###########
# functions
###########
# simBYM: simulates counts from BYM (IAR) model
# simLeroux: simulates counts from thinned Leroux model
# simUniform: simulate counts from uniform surveys


#########
## simBYM
#########
# simulates counts from BYM (IAR) model
# uses standardized (mean 0, variance 1) covariates
# returns a list of 
  # y: response count,
  # x: design matrix, intercept and covariate value
  # x.standardised: standardised x, x.sd: sd of covariate 
  # p: number of parameters 

## arguments
  # r: RasterBrick of covariates
  # beta: parameters vector, including intercept
  # sigma.sq: independent random effect variance
  # tau.sq: variance of the spatial random effects
  # rho: spatial strength of association parameter
  # nNeighbors: vector of # neighbors for each unit
simBYM <- function(r, beta, tau.sq, sigma.sq, nNeighbors){
  
  
  # spatial random effects
  D <- diag(nNeighbors)
  k <- nrow(D)
  Qw <- D - W
  Sigma <- tau.sq * ginv(Qw)
  phi <- mvrnorm(n=1, mu=rep(0, k), Sigma)
  
  
  # independent random effects
  theta <- rnorm(n=k, mean=0, sd=sqrt(sigma.sq))
  
  
  # covariates
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.i.standard <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }

  
  # simulate counts
  rates <- exp(x.standardised %*% beta + phi + theta)
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  sim <- list()
  sim$x <- x
  sim$x.standardised <- x.standardised
  sim$y <- counts
  sim$p <- length(beta)
  return(sim)
  
  
}


#########
## simGLM
#########
# simulates counts from Poisson glm
# uses standardized (mean 0, variance 1) covariates
# returns a list of y: response count,
  # x: design matrix, intercept and covariate value
  # x.standard: standardized x, x.sd: sd of covariate 
  # x.mean: mean of covariate
  # p: number  x parameters 

## arguments
  # r: RasterStack of covariates
  # beta: vector of parameters
simGLM <- function(r, beta){
  
  
  # covariates
  k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.i.standard <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }
  
  
  # simulate counts
  rates <- exp(x.standardised %*% beta)
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  sim <- list()
  sim$x <- x
  sim$x.standardised <- x.standardised
  sim$y <- counts
  sim$x.sd <- sd(x[,2])
  sim$x.mean <- mean(x[,2])
  sim$p <- 2
  return(sim)
  
  
}


############
## simLeroux
############
# simulates counts from Leroux model
# uses standardized (mean 0, variance 1) covariates
# returns a list of 
# y: response count,
# x: design matrix, intercept and covariate value
# x.standardised: standardised x, x.sd: sd of covariate 
# p: number of parameters

## params
# r: RasterBrick of covariates
# beta: parameters vector, including intercept
# tau.sq: variance of the spatial random effects
# rho: spatial strength of association parameter
# nNeighbors: vector of # neighbors for each unit
simLeroux <- function(r, beta, tau.sq, rho, nNeighbors, seed=NULL){
  
  
  # spatial random effects
  D <- diag(nNeighbors)
  k <- nrow(D)
  Qwp <- rho*(D - W) + (1-rho)*diag(k)
  Sigma <- tau.sq * ginv(Qwp)
  if (!is.null(seed)){ set.seed(seed) }
  phi <- mvrnorm(n=1, mu=rep(0, k), Sigma)
  
  
  # covariates
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.i.standard <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }
  
  
  # simulate counts
  rates <- exp(x.standardised %*% beta + phi)
  if (!is.null(seed)){ set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  sim <- list()
  sim$x <- x
  sim$x.standardised <- x.standardised
  sim$y <- counts
  sim$p <- length(beta)
  return(sim)
  
  
}


################
## simLerouxBias
################
# simulates counts from Leroux model
# uses standardized (mean 0, variance 1) covariates
# returns a list of 
# y: response count,
# x: design matrix, intercept and covariate value
# x.standardised: standardised x, x.sd: sd of covariate 
# p: number of parameters

## params
# r: RasterBrick of covariates
# beta: parameters vector, including intercept
# gamma: bias (thinning) parameter
# tau.sq: variance of the spatial random effects
# rho: spatial strength of association parameter
# nNeighbors: vector of # neighbors for each unit
simLerouxBias <- function(r, beta, gamma, tau.sq, rho, nNeighbors, seed=NULL){
  
  
  # spatial random effects
  D <- diag(nNeighbors)
  k <- nrow(D)
  Qwp <- rho*(D - W) + (1-rho)*diag(k)
  Sigma <- tau.sq * ginv(Qwp)
  if (!is.null(seed)){ set.seed(seed) }
  phi <- mvrnorm(n=1, mu=rep(0, k), Sigma)
  
  
  # covariates
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.i.standard <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }
  
  
  # simulate counts
  rates <- exp(x.standardised %*% beta + gamma + phi)
  if (!is.null(seed)){ set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  sim <- list()
  sim$x <- x
  sim$x.standardised <- x.standardised
  sim$y <- counts
  sim$p <- length(beta)
  return(sim)
  
  
}


####################
## simLerouxThinned
###################
# simulates counts from thinned Leroux model
# uses standardized (mean 0, variance 1) covariates
# returns a list of 
  # y: response count,
  # x: design matrix, intercept and covariate value
  # x.standardised: standardised x, x.sd: sd of covariate 
  # p: number of parameters

## arguments
  # r: RasterBrick of covariates
  # r.thin: RasterBrick of thinning covariates
  # beta: covariate params, including intercept
  # gamma: vector of thinning covariate parameters
  # tau.sq: variance of the spatial random effects
  # rho: spatial strength of association parameter
  # nNeighbors: vector of # neighbors for each unit
simLerouxThinned <- function(r, r.thin, beta, gamma, tau.sq, rho, nNeighbors, seed=NULL){
  
  
  # spatial random effects
  D <- diag(nNeighbors)
  k <- nrow(D)
  Qwp <- rho*(D - W) + (1-rho)*diag(k)
  Sigma <- tau.sq * ginv(Qwp)
  if (!is.null(seed)){ set.seed(seed) }
  phi <- mvrnorm(n=1, mu=rep(0, k), Sigma)

  
  # covariates
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.i.standard <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }
  
  
  # thinning covariates
  z <- array(NA, c(k, length(names(r.thin))))
  z.standardised <- z
  for (i in 1:length(names(r.thin))){
    z.i <- values(r.thin[[i]])
    z.i <- z.i[!is.na(z.i)]
    z.i.standard <- (z.i - mean(z.i))/sd(z.i)
    z[,i] <- z.i
    z.standardised[,i] <- z.i.standard
  }
  
  
  # simulate counts
  rates <- exp(x.standardised %*% beta + z.standardised %*% gamma + phi)
  counts <- sapply(rates, function(x){ 
    if (!is.null(seed)){ set.seed(seed) }
  rpois(n=1, x)})
  
  
  # create the response
  sim <- list()
  sim$x <- x
  sim$x.standardised <- x.standardised
  sim$z <- z
  sim$z.standardised <- z.standardised
  sim$y <- counts
  sim$p <- length(beta)
  sim$p.z <- length(gamma)
  return(sim)
  
  
}


simUniform <- function(r, beta, nsamp, seed=NULL){
  
  
  # covariates
  k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.i.standard <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }
  
  
  # simulate counts
  rates <- exp(x.standardised %*% beta)
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # sample locations
  if (!is.null(seed)){ set.seed(seed) }
  cell.ids <- sample(k, nsamp)
  samp.counts <- counts[cell.ids]
  samp.x <- x[cell.ids,]
  samp.x.stand <- x.standardised[cell.ids,]
  
  
  # create the response
  sim <- list()
  sim$x <- samp.x
  sim$x.standardised <- samp.x.stand
  sim$y <- samp.counts
  sim$x.sd <- sd(x[,2])
  sim$x.mean <- mean(x[,2])
  sim$p <- length(beta)
  sim$cell.ids <- cell.ids
  return(sim)
  
  
}


## simBernoulliLoc: simulates survey locations
  # calculates survey occupancy status as bernoulli RV
  # whose success probability is a function of a 
  # covariate surface.
# returns a vector of occupancy indicators 
  # corresponding to raster cell ids

## arguments
  # r: raster of covariates governing sampling location
  # beta: sampling parameters
simBernoulliLoc <- function(r, beta, seed=NULL){
  
  
  # covariates
  if (length(beta) == length(names(r))){
    
    # no intercept, scaled
    k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
    x <- array(NA, c(k, length(names(r))))
    x.scaled <- x
    for (i in 1:length(names(r))){
      x.i <- values(r[[i]])
      x.i <- x.i[!is.na(x.i)]
      x.scaled.i <- x.i/sd(x.i)
      x[,i] <- x.i
      x.scaled[,i] <- x.scaled.i
    }
    
  } else {
    
    # intercept, standardised
    k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
    x <- array(1, c(k, 1))
    x.scaled <- x
    for (i in 1:length(names(r))){
      x.i <- values(r[[i]])
      x.i <- x.i[!is.na(x.i)]
      x.scaled.i <- (x.i - mean(x.i))/sd(x.i)
      x <- cbind(x, x.i)
      x.scaled <- cbind(x.scaled, x.scaled.i)
    }
    
  }
  
  
  # simulate survey occupancy
  lin.pred <- x.scaled %*% beta
  probs <- exp(lin.pred)/(1 + exp(lin.pred))
  if (!is.null(seed)){ set.seed(seed) }
  status <- sapply(probs, function(x){
    rbinom(n=1, size=1, x)})
  
  
  cells.all <- c(1:ncell(r))[!is.na(values(r[[1]]))]
  cells <- cells.all[as.logical(status)]
  output <- list()
  output$status <- status
  output$cells <- cells
  output$coords <- xyFromCell(loc.disc, cell=cells)
  output$probs <- probs
  output$p <- ncol(x)
  output$x <- x
  output$x.scaled <- x.scaled
  return(output)
  
}


# simulates survey locations according 
# to covariates and a Gaussian process
simBernoulliLocGp <- function(r, beta, seed=NULL, w=NULL){
  
  
  # covariates
  if (length(beta) == length(names(r))){
    
    # no intercept, scaled
    k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
    x <- array(NA, c(k, length(names(r))))
    x.scaled <- x
    for (i in 1:length(names(r))){
      x.i <- values(r[[i]])
      x.i <- x.i[!is.na(x.i)]
      x.scaled.i <- x.i/sd(x.i)
      x[,i] <- x.i
      x.scaled[,i] <- x.scaled.i
    }
    
  } else {
    
    # intercept, standardised
    k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
    x <- array(1, c(k, 1))
    x.scaled <- x
    for (i in 1:length(names(r))){
      x.i <- values(r[[i]])
      x.i <- x.i[!is.na(x.i)]
      x.scaled.i <- (x.i - mean(x.i))/sd(x.i)
      x <- cbind(x, x.i)
      x.scaled <- cbind(x.scaled, x.scaled.i)
    }
    
  }
  
  
  # simulate survey occupancy
  lin.pred <- x.scaled %*% beta + w
  probs <- exp(lin.pred)/(1 + exp(lin.pred))
  if (!is.null(seed)){ set.seed(seed) }
  status <- sapply(probs, function(x){
    rbinom(n=1, size=1, x)})
  
  
  # ids of observed cells
  ids <- c()
  for (i in 1:length(status)){
    if (as.logical(status[i])){
      ids <- c(ids, i)
    }
  }
  
  
  cells <- cells.all[as.logical(status)]
  output <- list()
  output$status <- status
  output$cells <- cells
  output$coords <- xyFromCell(loc.disc, cell=cells)
  output$probs <- probs
  output$p <- ncol(x)
  output$x <- x
  output$x.scaled <- x.scaled
  output$w <- w
  output$d <- d
  output$ids <- ids
  return(output)
  
  
}


simBernoulliLocCov <- function(r, beta, w, seed=NULL){
  
  
  # intercept, standardised
  k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
  x <- array(1, c(k, 1))
  x.scaled <- x
  for (i in 1:length(names(r))){
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.scaled.i <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.scaled <- cbind(x.scaled, x.scaled.i)
  }
    
  
  # simulate survey occupancy
  lin.pred <- x.scaled %*% beta + w
  probs <- exp(lin.pred)/(1 + exp(lin.pred))
  if (!is.null(seed)){ set.seed(seed) }
  status <- sapply(probs, function(x){
    rbinom(n=1, size=1, x)})
  
  
  # ids of observed cells
  ids <- c()
  for (i in 1:length(status)){
    if (as.logical(status[i])){
      ids <- c(ids, i)
    }
  }
  
  
  cells <- cells.all[as.logical(status)]
  output <- list()
  output$status <- status
  output$cells <- cells
  output$coords <- xyFromCell(r[[1]], cell=cells)
  output$probs <- probs
  output$p <- ncol(x)
  output$x <- x
  output$x.scaled <- x.scaled
  output$w <- w
  output$d <- d
  output$ids <- ids
  return(output)
  
  
}


# to covariates and a Gaussian process
simLocW <- function(w, r, beta=0, seed=NULL){
  
  
  # simulate survey occupancy
  lin.pred <- beta + w
  probs <- exp(lin.pred)/(1 + exp(lin.pred))
  if (!is.null(seed)){ set.seed(seed) }
  status <- sapply(probs, function(x){
    rbinom(n=1, size=1, x)})
  
  
  # ids of observed cells
  ids <- c()
  for (i in 1:length(status)){
    if (as.logical(status[i])){
      ids <- c(ids, i)
    }
  }
  
  
  cells <- cells.all[as.logical(status)]
  output <- list()
  output$status <- status
  output$cells <- cells
  output$coords <- xyFromCell(r, cell=cells)
  output$probs <- probs
  output$ids <- ids
  return(output)
  
  
  
}


simConditional <- function(r, loc.stats, beta, beta.samp, alpha, global.center=FALSE, seed=NULL){
  
  
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
  
  
  # extract cell selection probability
  probs <- loc.stats$probs[as.logical(loc.stats$status)]

  
  # simulate counts
  rates <- exp(x.standardised %*% beta + alpha * probs)
  if (!is.null(seed)){ set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  output <- list()
  output$y <- counts
  output$x.standardised <- x.standardised
  output$x <- x
  output$p <- ncol(x)
  return(output)
    
  
}


simConditionalGp <- function(r, loc.stats, beta, beta.samp, alpha, global.center=FALSE, seed=NULL){
  
  
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
  
  
  # extract cell selection probability
  probs <- loc.stats$probs[as.logical(loc.stats$status)]
  
  
  # simulate counts
  rates <- exp(x.standardised %*% beta + alpha * probs)
  if (!is.null(seed)){ set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  output <- list()
  output$y <- counts
  output$x.standardised <- x.standardised
  output$x <- x
  output$p <- ncol(x)
  return(output)
  
  
}


simConditionalGp2 <- function(r, loc.stats, beta, alpha, w, center=TRUE, global.center=FALSE, seed=NULL, offset=1){
  
  
  # extract covariate values at survey locations
  k <- length(loc.stats$cells)
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- raster::extract(r[[i]], loc.stats$coords)
    if (global.center){
      vals <- values(r[[i]])[!is.na(values(r[[i]]))]
      mu <- mean(vals)
      sdev <- sd(vals) 
    } else if (center){
      mu <- mean(x.i)
      sdev <- sd(x.i)
    } else{
      mu <- 0
      sdev <- 1
    }
    x.i.standard <- (x.i - mu)/sdev
    x <- cbind(x, x.i)
    x.standardised <- cbind(x.standardised, x.i.standard)
  }
  
  
  # extract cell selection probability
  w.sub <- w[as.logical(loc.stats$status)]
  
  
  # simulate counts
  rates <- exp(log(offset) + x.standardised %*% beta + alpha * w.sub)
  if (!is.null(seed)){ set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  output <- list()
  output$y <- counts
  output$x.standardised <- x.standardised
  output$x <- x
  output$p <- ncol(x)
  return(output)
  
  
}


simConditionalMVGP <- function(r, loc.stats, beta, alpha, w, g, global.center=FALSE, seed=NULL, offset=1){
  
  
  # extract covariate values at survey locations
  k <- length(loc.stats$cells)
  x <- matrix(rep(1, k), ncol=1)
  x.standardised <- x
  for (i in 1:length(names(r))){
    x.i <- raster::extract(r[[i]], loc.stats$coords)
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
  
  
  # extract cell selection probability
  w.sub <- w[as.logical(loc.stats$status)]
  g.sub <- g[as.logical(loc.stats$status)]
  
  # simulate counts
  rates <- exp(log(offset) + x.standardised %*% beta + alpha * w.sub + g.sub)
  if (!is.null(seed)){ set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n=1, x)})
  
  
  # create the response
  output <- list()
  output$y <- counts
  output$x.standardised <- x.standardised
  output$x <- x
  output$p <- ncol(x)
  return(output)
  
  
}


#############
## locCondSim
#############
# simulates location conditional data: 
# i.e. survey locations, and disease
# counts given location. 
# returns a list with elements:
# conditional: output of simConditional
# loc: output of simBernoulliLoc

## arguments
# r: raster of disease related covariates
# r.samp: raster of sampling location covariates
# beta: parameter for disease covariates
# beta.samp: parameter for sampliong covariates
# alpha: sampling scale parameter
simLocCond <- function(r, r.samp, beta, beta.samp, alpha, global.center=FALSE){
  
  
  #### Simulate survey locations
  locs <- simBernoulliLoc(r.samp, beta.samp)
  
  
  #### Simulate counts given locations
  count.data <- simConditional(r, locs, beta, beta.samp, alpha, global.center)
  data <- list(conditional=count.data, loc=locs)
  return(data)
  
  
}


# simulate CAR random effects
simLerouxThinned <- function(tau.sq, rho, nNeighbors, seed=NULL){
  
  
  D <- diag(nNeighbors)
  k <- nrow(D)
  Qwp <- rho*(D - W) + (1-rho)*diag(k)
  Sigma <- tau.sq * ginv(Qwp)
  if (!is.null(seed)){ set.seed(seed) }
  phi <- mvrnorm(n=1, mu=rep(0, k), Sigma)
  return(phi)
  

}
