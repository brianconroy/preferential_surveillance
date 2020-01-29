library(fields)
library(MASS)


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
