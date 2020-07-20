library(fields)
library(MASS)
library(raster)

simLocations <- function(r, beta, w, cells.all, d, seed=NULL){
  
  
  # intercept, standardized
  k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
  x <- array(1, c(k, 1))
  x.scaled <- x
  for (i in 1:length(names(r))) {
    x.i <- values(r[[i]])
    x.i <- x.i[!is.na(x.i)]
    x.scaled.i <- (x.i - mean(x.i))/sd(x.i)
    x <- cbind(x, x.i)
    x.scaled <- cbind(x.scaled, x.scaled.i)
  }
    
  
  # simulate survey occupancy
  lin.pred <- x.scaled %*% beta + w
  probs <- exp(lin.pred)/(1 + exp(lin.pred))
  if (!is.null(seed)) { set.seed(seed) }
  status <- sapply(probs, function(x){
    rbinom(n = 1, size = 1, x)})
  
  
  # ids of observed cells
  ids <- c()
  for (i in 1:length(status)) {
    if (as.logical(status[i])) {
      ids <- c(ids, i)
    }
  }
  
  
  cells <- cells.all[as.logical(status)]
  output <- list()
  output$status <- status
  output$cells <- cells
  output$coords <- xyFromCell(r[[1]], cell = cells)
  output$probs <- probs
  output$p <- ncol(x)
  output$x <- x
  output$x.scaled <- x.scaled
  output$w <- w
  output$d <- d
  output$ids <- ids
  return(output)
  
  
}
