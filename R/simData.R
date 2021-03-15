

simDataset <- function(r, beta.loc, beta.case, beta.ctrl, alpha.case, alpha.ctrl,
                       theta=NULL, phi=NULL, w=NULL){
  
  
  # get ids of non-null cells in input raster
  cells.all <- c(1:ncell(r))[!is.na(values(r[[1]]))]
  d <- form_distance_matrix(r, cells.all)
  
  if (is.null(w)){
    # simulate gaussian process
    sigma <- Exponential(d, range = theta, phi = phi)
    w <- mvrnorm(n = 1, mu = rep(0, length(cells.all)), sigma)
  }
  
  # simulate locations
  locs <- simLocations(r, beta.loc, w = w, cells.all = cells.all, d = d)
  
  # Simulate counts given locations
  case.data <- simCounts(r, locs, beta.case, alpha.case, w)
  ctrl.data <- simCounts(r, locs, beta.ctrl, alpha.ctrl, w)
  
  return(list(
    locs=locs,
    case.data=case.data,
    ctrl.data=ctrl.data
  ))
  
  
}


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


#' simLocations
#' 
#' In this package, a disease surveillance system collects samples for a disease
#' over cells in a raster. The sampling process may be preferential,
#' assigning sampling locations in a way that is stochastically related
#' to the disease process. Typically preferential sampling manifests itself when
#' sampling locations are generally assigned in areas thought to be at high risk
#' for the disease.
#' 
#' simLocations simulates sample sites preferentially over a raster study region.
#'
#' @param r: RasterBrick. Raster representing the study region whose values are covariates associated with 
#' the sampling process.
#' @param beta: numeric. Vector of parameters (including intercept) associated with the covariates in r. 
#' @param w: numeric. Vector of spatial random effects simulated from a Gaussian process.  
#' @param cells.all: numeric. Vector of integers representing id's of non-null cells in raster r.
#' @param d: matrix. Distance matrix describing distances between cells in raster r. 
#' @param seed: numeric. Optional random seed value for reproducibility.
#'
#' @return a list with the following keys:
#' 
#' \itemize{
#'   \item status: numeric. Vector of 0/1 indicators representing whether the corresponding cell 
#'   in cells.all was observed.
#'   \item cells: numeric. Vector of cell ids from cells.all that were observed.
#'   \item coords: matrix. Matrix of latitude (x) and longitude (y) coordinates from observed cells.
#'   \item x: matrix. Matrix of intercept and covariate values from raster r at the center point of 
#'   each cell in r.
#' }
#' 
#' @export
#'
#' @examples
#' 
#' # create covariate raster
#' r <- aggregate(load_prism_pcs(), fact=10)
#' 
#' # get cell ids
#' cells.all <- c(1:ncell(r))[!is.na(values(r[[1]]))]
#' 
#' # get distance matrix 
#' d <- form_distance_matrix(r, cells.all)
#' 
#' # simulate Gaussian process
#' sigma <- Exponential(d, range = 7, phi = 5)
#' w <- mvrnorm(n = 1, mu = rep(0, length(cells.all)), sigma)
#' 
#' # simulate locations
#' locs <- simLocations(r, c(-1.50, 1.00, -0.25), w, cells.all, d)
#' 
#' # view raster of observed cells
#' r_loc <- r[[1]]
#' r_loc[!is.na(r_loc[])] <- locs$status
#' plot(r_loc)
#' 
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


form_distance_matrix <- function(r, cells.all){
  
  coords <- xyFromCell(r, cell=cells.all)
  return(as.matrix(dist(coords, diag=TRUE, upper=TRUE)))
  
}
