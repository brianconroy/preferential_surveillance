library(fields)
library(MASS)
library(raster)


simCounts <- function(r, loc.stats, beta, alpha, w, center=TRUE, global.center=FALSE, seed=NULL, offset=1){
  
  
  # extract covariate values at survey locations
  k <- length(loc.stats$cells)
  x <- matrix(rep(1, k), ncol = 1)
  x.standardised <- x
  for (i in 1:length(names(r))) {
    x.i <- raster::extract(r[[i]], loc.stats$coords)
    if (global.center) {
      vals <- values(r[[i]])[!is.na(values(r[[i]]))]
      mu <- mean(vals)
      sdev <- sd(vals) 
    } else if (center) {
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
  if (!is.null(seed)) { set.seed(seed) }
  counts <- sapply(rates, function(x){rpois(n = 1, x)})
  
  
  # create the response
  output <- list()
  output$y <- counts
  output$x.standardised <- x.standardised
  output$x <- x
  output$p <- ncol(x)
  return(output)
}
