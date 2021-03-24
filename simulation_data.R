###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (spatial poisson regression,
# poisson regression, and
# BART classifiers) under
# low and high levels of 
# preferential sampling.

# This script simulated datasets
# under low and high levels of
# preferential sampling.
###############################

R.utils::sourceDirectory('R/')

set.seed(12)

# Prism principal components ----------------------------------------------

agg_factor <- 10
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact = agg_factor)
print(n_values(caPr.disc[[1]]))
plot(caPr.disc)

# Simulate data function --------------------------------------------------

simulate_data <- function(label, disc, n_sims=25, save_dir=NULL, 
                          Alpha.case=NULL, Alpha.ctrl=NULL, 
                          beta.case=NULL, beta.ctrl=NULL, beta.loc=NULL,
                          Theta=6, Phi=12){
  
  # Alpha.case: represents the strength of preferential sampling on 
  # observed distributions of cases (disease positive specimen). 
  # The greater Alpha.case, the stronger the dependence between 
  # sampling location and the abundance of cases. 
  
  # Alpha.ctrl: same interpretation as Alpha.case but for controls
  # (disease negative specimen).

  # beta.case: parameters for fixed effect covariates associated w/ the spatial
  # distribution of cases. The vector consists of an intercept 
  # (first element) and as many additional parameters as there are 
  # layers in the "disc" raster
  
  # beta.ctrl: same as above but for controls.
  
  # beta.loc: parameters for fixed effect covariates associated with the spatial
  # distribution of observation sites.
  
  if (label == "low") {
    Alpha.case <- 0.5 
    Alpha.ctrl <- 0.3
    beta.case <- c(1, 0.75, 0.25)
    beta.ctrl <- c(3, 1, 0.5)
    beta.loc <- c(-1.5, 1, -0.25)
  } else if (label == "high") {
    Alpha.case <- 1
    Alpha.ctrl <- -0.25
    beta.case <- c(-1.5, 0.25, 0.25)
    beta.ctrl <- c(3, 1, 0.5)
    beta.loc <- c(-1.5, 1, -0.25)
  }
  
  prevalences <- c()
  ps_contribs <- c()
  obs_cells <- c()
  n_specimen <- c()
  
  for (i in 1:n_sims) {
    
    #### Simulate Gaussian process
    cells.all <- c(1:ncell(disc))[!is.na(values(disc[[1]]))]
    coords <- xyFromCell(disc, cell = cells.all)
    d <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
    Sigma <- Exponential(d, range = Theta, phi = Phi)
    W <- mvrnorm(n = 1, mu = rep(0, length(cells.all)), Sigma)
    
    #### Simulate locations
    locs <- simLocations(disc, beta.loc, w = W, cells.all = cells.all, d = d)
    obs_cells <- c(obs_cells, sum(locs$status))
    
    # the measure of the effects of preferential sampling on observed risk
    ps_contribs <- c(ps_contribs, 
                     calc_ps_contribution(disc, locs, beta.case, Alpha.case, beta.ctrl, Alpha.ctrl, W))
    
    #### Simulate counts given locations
    case.data <- simCounts(disc, locs, beta.case, Alpha.case, W)
    ctrl.data <- simCounts(disc, locs, beta.ctrl, Alpha.ctrl, W)
    prevalences <- c(prevalences, sum(case.data$y)/sum(case.data$y + ctrl.data$y))
    n_specimen <- c(n_specimen, sum(case.data$y + ctrl.data$y))
    
    params <- list(
      sampling = label,
      beta.case = beta.case,
      beta.ctrl = beta.ctrl,
      alpha.case = Alpha.case,
      alpha.ctrl = Alpha.ctrl,
      beta.loc = beta.loc,
      Theta = Theta,
      Phi = Phi,
      W = W
    )
    
    data <- list(
      case.data = case.data,
      ctrl.data = ctrl.data,
      locs = locs
    )
    
    if (!is.null(save_dir)) {
      dir.create(save_dir)
      param_fname <- paste0("params_", label, "_", i, ".json")
      data_fname  <- paste0("data_", label, "_", i, ".json")
      
      write( toJSON(params), file.path(save_dir, param_fname) )
      write( toJSON(data), file.path(save_dir, data_fname) )
    }
  }
  
  list(prevalences = prevalences, 
       ps_contribs = ps_contribs, 
       obs_cells = obs_cells,
       n_specimen = n_specimen)
}

# Generate high and low preferential samples ------------------------------

low_sample <- simulate_data("low", caPr.disc)
high_sample <- simulate_data("high", caPr.disc)
