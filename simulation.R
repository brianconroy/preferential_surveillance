###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (spatial poisson regression,
# poisson regression, and
# BART classifiers) under
# low and high levels of 
# preferential sampling.

# This script fits models to
# the already simulated datasets.
###############################


library(raster)
library(sp)
library(plyr)
library(mvtnorm)
library(MASS)
library(fields)
library(R.utils)
library(dbarts)
library(preferentialSurveillance)
#sourceDirectory('R/')


#############################################################################
#                       Low Preferential Sampling                           #
#############################################################################

# The first set of simulations is run under a "low" level of preferential
# sampling, i.e. a generally weaker association between the location of 
# sample sites and underlying disease risk.

# Setup ---------------------------------------------------------------------

sampling <- "low"
src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")
n_sims <- 25
agg_factor <- 10

#### Prism Principal Components
caPr <- prism_pca
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))

###################################################
#                Proposed Method                  #
###################################################

for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  
  # Estimate initial values for spatial params ------------------------------
  params <- load_output(paste("params_low_", i, ".json", sep=""), src=src)
  
  # Specify priors for fitting initial values. 
  # These priors are NOT used when fitting 
  # the full model.
  Theta <- params$Theta
  Phi <- params$Phi
  prior_theta <- get_gamma_prior(Theta, 5)
  prior_phi <- get_igamma_prior(Phi, 5)
  
  # Estimate an initial value for the vector
  # of spatial random effects (W). Do so by
  # regressing binary indicators for whether
  # a grid cell was observed on prism covariates
  # plus realizations from a Gaussian process (W)
  w_output <- logisticGp(
                         # binary indicators, covariates, distance matrix
                         y=data$locs$status, x=data$locs$x.scaled, d,
                         
                         # MCMC sampling params
                         n.sample=1500, burnin=500, L_beta=8, L_w=8, proposal.sd.theta=0.3,
                         
                         # Initial values
                         w_initial=NULL, theta_initial=NULL, phi_initial=NULL, beta_initial=NULL,
                         
                         # Prior parametrizations
                         prior_phi=prior_phi, prior_theta=prior_theta)
  
  # Specify initial MCMC values
  w_initial <- colMeans(w_output$samples.w)
  theta_initial <- mean(w_output$samples.theta)
  phi_initial <- mean(w_output$samples.phi)
  beta_loc_initial <- colMeans(w_output$samples.beta)
  
  # Estimate beta+ and alpha+ initial values ------------------------------
  ini_case <- glm(data$case.data$y ~ data$case.data$x.standardised + w_initial[data$locs$ids] - 1, family='poisson')
  alpha_ca_initial <- coefficients(ini_case)[4]
  beta_ca_initial <- coefficients(ini_case)[1:3]
  
  # Estimate beta- and alpha- initial values ------------------------------
  ini_ctrl <- glm(data$ctrl.data$y ~ data$ctrl.data$x.standardised + w_initial[data$locs$ids] - 1, family='poisson')
  alpha_co_initial <- coefficients(ini_ctrl)[4]
  beta_co_initial <- coefficients(ini_ctrl)[1:3]
  
  # Specify priors --------------------------------------------------------
  prior_alpha_ca_mean <- alpha_ca_initial
  prior_alpha_ca_var <- 3
  prior_alpha_co_mean <- alpha_co_initial
  prior_alpha_co_var <- 3
  prior_theta <- get_gamma_prior(theta_initial, 2)
  prior_phi <- get_igamma_prior(phi_initial, 2)
  
  # MCMC params ------------------------------------------------------------
  n.sample <- 10000
  burnin <- 3000
  proposal.sd.theta <- 0.15
  
  # Hamiltonian Monte Carlo parameters
  L_w <- 8
  L_ca <- 8
  L_co <- 8
  L_a_ca <- 8
  L_a_co <- 8
  L_loc <- 8
  
  # Dual averaging parameters
  m_aca <- 2000
  m_aco <- 2000
  m_ca <- 2000
  m_co <- 2000
  m_w <- 2000
  m_loc <- 2000
  target_aca=0.65
  target_aco=0.65
  target_ca=0.65
  target_co=0.65
  target_w=0.65
  target_loc=0.65
  
  # Run fit
  output <- preferentialSampling(
                                # data, distance matrix
                                data, d, 
                                
                                # MCMC params
                                n.sample, burnin,
                                L_w, L_ca, L_co, L_a_ca, L_a_co,
                                proposal.sd.theta=proposal.sd.theta,
                           
                                # HMC self tuning params
                                m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, 
                                target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, 
                                target_co=target_co, target_w=target_w, target_loc=target_loc,
                                self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, 
                                self_tune_co=TRUE, self_tune_loc=TRUE,
                                
                                # Initial MCMC values
                                beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, 
                                alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial, 
                                beta_loc_initial=beta_loc_initial, theta_initial=theta_initial, 
                                phi_initial=phi_initial, w_initial=w_initial,
                                
                                # Priors
                                prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, 
                                prior_alpha_co_var=prior_alpha_co_var
                                )
  
  # Check estimated log odds -------------------------------------------------
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  # True param values
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  # True log odds vs estimated
  X.standard <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  
  # Save output --------------------------------------------------------------
  output$description <- paste(sim_name, "_", i, sep="")
  save_output(output, paste("output_", sim_name, "_", i, ".json", sep=""), dst=src)
  
}


###################################################
#           Spatial Poisson Regression            #
###################################################

for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_low_", i, ".json", sep=""), src=src)
  
  #########
  ### Cases
  #########
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  d.sub <- d[as.logical(data$locs$status), as.logical(data$locs$status)]
  
  ## Fit model for cases
  # Initial values
  beta_ca_i_ <- rnorm(ncol(X.ca))
  w_i_ <- rnorm(nrow(d.sub))
  phi_i_ <- params$Phi + rnorm(1)
  theta_i_ <- params$Theta + rnorm(1)
  
  # Other parameters
  n.sample_ <- 10000
  burnin_ <- 3000
  L_w_ <- 8
  L_b_ <- 8
  prior_phi_ <- c(3, 40)
  prior_theta_ <- get_gamma_prior(prior_mean=params$Theta, prior_var=5)
  
  # Run fit
  output.sp_ca <- poissonGp(X.ca, Y.ca, d.sub,
                            n.sample=n.sample_, burnin=burnin_, proposal.sd.theta=0.3,
                            L_w=L_w_, L_b=L_b_,
                            beta_initial=beta_ca_i_, w_initial=w_i_, 
                            phi_initial=phi_i_, theta_initial=theta_i_,
                            prior_phi=prior_phi_, prior_theta=prior_theta_)
  
  # Krige random effects
  w.hat_spca <- colMeans(output.sp_ca$samples.w)
  
  kriged_w_ca <- krigeW(output.sp_ca, d, data$locs$ids)
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(data$locs$status))
  
  # Save output and kriged estimates
  output.sp_ca$description <- paste(sim_name, "_", i, sep="")
  save_output(output.sp_ca, paste("output.sp_ca_", sim_name, "_", i, ".json", sep=""), dst=src)
  save_output(kriged_w_ca, paste("output.krige_ca_", sim_name, "_", i, ".json", sep=""), dst=src)
  
  ############
  ### Controls
  ############
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  
  n.sample__ <- 10000
  burnin__ <- 3000
  L_w__ <- 8
  L_b__ <- 8
  
  beta_co_i__ <- rnorm(ncol(X.co))
  w_i__ <- rnorm(nrow(d.sub))
  phi_i__ <- params$Phi + rnorm(1)
  theta_i__ <- params$Theta + rnorm(1)
  
  prior_phi__ <- c(3, 40)
  prior_theta__ <- c(2.5, 2.5)
  
  output.sp_co <- poissonGp(X.co, Y.co, d.sub,
                            n.sample=n.sample__, burnin=burnin__, 
                            L_w=L_w__, L_b=L_b__, proposal.sd.theta=0.3,
                            beta_initial=beta_co_i__, w_initial=w_i__, 
                            phi_initial=phi_i__, theta_initial=theta_i__,
                            prior_phi=prior_phi__, prior_theta=prior_theta__)
  
  # Krige random effects
  w.hat_spco <- colMeans(output.sp_co$samples.w)
  
  kriged_w_co <- krigeW(output.sp_co, d, data$locs$ids)
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(data$locs$status))
  
  output.sp_co$description <- paste(sim_name, "_", i, sep="")
  save_output(output.sp_co, paste("output.sp_co_", sim_name, "_", i, ".json", sep=""), dst=src)
  save_output(kriged_w_co, paste("output.krige_co_", sim_name, "_", i, ".json", sep=""), dst=src)
  
}


######################################
# Bayesian Additive Regression Trees #
# (Classifier)                       #
######################################

# Probit BART classifiers are fit to the (transformed) simulated
# datasets. We transform our simulated datasets of aggregated counts 
# into datasets of binary positive-negative outcomes. Similarly to 
# Carlson et al. (2021), we use the observed Y_i+ positive counts 
# in a grid cell to generate a set of Y_i+ pseudo-positive points 
# located at the centroid of the cell, and Y_i- counts of negative
# specimen to generate a set of Y_i- pseudo-negative points in the 
# same location

for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  
  #########
  ### Cases
  #########
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  X.coords <- data$locs$coords
  X.combined <- cbind(X.ca, X.coords)
  
  # Generate pseudo-positive points
  Y.new <- c()
  X.new.combined <- c()
  for (idx in 1:length(Y.ca)){
    
    num_ca <- Y.ca[idx]
    Y.new <- c(Y.new, rep(1, num_ca))
    
    # Repeat the covariates for the grid cell for each point
    if (num_ca > 0){
      for (iter in 1:num_ca){
        X.new.combined <- rbind(X.new.combined, X.combined[idx,])
      }
    }
  
  }
  
  ############
  ### Controls
  ############
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  
  for (idx in 1:length(Y.co)){
    
    num_co <- Y.co[idx]
    Y.new <- c(Y.new, rep(0, num_co))
    
    if (num_co > 0){
      for (iter in 1:num_co){
        X.new.combined <- rbind(X.new.combined, X.combined[idx,])
      }
    }
    
  }
  
  #######################
  ### Fit Bart Classifier
  #######################
  
  # By default a probit link function is used if Y has only values 0, 1
  bart_clsf <- bart(X.new.combined, Y.new,
                  # default prior hyperparams from the original BART paper
                  k = 2.0,
                  power = 2.0, 
                  base = 0.95,
                  # as suggested by BART authors
                  ntree = 200,
                  # number of posterior samples
                  ndpost = 2000, 
                  # burnin
                  nskip = 500,
                  printevery = 100, 
                  # no thinning
                  keepevery = 1, 
                  keeptrainfits = TRUE,
                  usequants = FALSE, 
                  numcut = 100, 
                  verbose = TRUE, 
                  # 4 independent MCMC chains
                  nchain = 4, 
                  nthread = 1, 
                  combinechains = TRUE,
                  keeptrees = TRUE)
  
  # Load covariates for all cells in study region
  X.full <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  
  # Combine with cordinates over the entire study region
  r <- caPr.disc[[1]]
  coords_all <- xyFromCell(caPr.disc[[1]], cell = (1:length(r[]))[!is.na(r[])])
  X.full.comb <- cbind(X.full, coords_all)
  
  # Estimate posterior over full study region
  posterior_risk <- predict(bart_clsf, X.full.comb)
  post_mean_risk <- colMeans(posterior_risk)
  
  # Histogram of risk
  hist(post_mean_risk, main = paste("Iteration", i))
  
  ## Save data ------------------------------------------------
  save_output(post_mean_risk, paste("output.bart_clsf_", sim_name, "_", i, ".json", sep=""), dst=src)

}

#############################################################################
#                       High Preferential Sampling                          #
#############################################################################

sampling <- "high"
src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")
n_sims <- 25
agg_factor <- 10

#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))

for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  
  ## Initial Values
  # Random effects, beta.samp, theta, phi
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
  Theta <- params$Theta
  Phi <- params$Phi
  prior_theta <- get_gamma_prior(Theta, 5)
  prior_phi <- get_igamma_prior(Phi, 5)
  w_output <- logisticGp(y=data$locs$status, x=data$locs$x.scaled, d, n.sample=1500, burnin=500, L_beta=8, L_w=8, proposal.sd.theta=0.3,
                         w_initial=NULL, theta_initial=NULL, phi_initial=NULL, beta_initial=NULL,
                         prior_phi=prior_phi, prior_theta=prior_theta)
  
  w_initial <- colMeans(w_output$samples.w)
  theta_initial <- mean(w_output$samples.theta)
  phi_initial <- mean(w_output$samples.phi)
  beta_loc_initial <- colMeans(w_output$samples.beta)
  
  # Beta & alpha (case) initial values
  ini_case <- glm(data$case.data$y ~ data$case.data$x.standardised + w_initial[data$locs$ids] - 1, family='poisson')
  alpha_ca_initial <- coefficients(ini_case)[4]
  beta_ca_initial <- coefficients(ini_case)[1:3]
  
  # Beta & alpha (control) initial values
  ini_ctrl <- glm(data$ctrl.data$y ~ data$ctrl.data$x.standardised + w_initial[data$locs$ids] - 1, family='poisson')
  alpha_co_initial <- coefficients(ini_ctrl)[4]
  beta_co_initial <- coefficients(ini_ctrl)[1:3]
  
  # Fit full model
  prior_alpha_ca_mean <- alpha_ca_initial
  prior_alpha_ca_var <- 3
  prior_alpha_co_mean <- alpha_co_initial
  prior_alpha_co_var <- 3
  prior_theta <- get_gamma_prior(theta_initial, 2)
  prior_phi <- get_igamma_prior(phi_initial, 2)
  
  n.sample <- 10000
  burnin <- 3000
  proposal.sd.theta <- 0.15
  
  # Hamiltonian Monte Carlo parameters
  L_w <- 8
  L_ca <- 8
  L_co <- 8
  L_a_ca <- 8
  L_a_co <- 8
  L_loc <- 8
  
  # Dual averaging parameters
  m_aca <- 2000
  m_aco <- 2000
  m_ca <- 2000
  m_co <- 2000
  m_w <- 2000
  m_loc <- 2000
  target_aca=0.65
  target_aco=0.65
  target_ca=0.65
  target_co=0.65
  target_w=0.65
  target_loc=0.65
  
  # Run fit
  output <- preferentialSampling(data, d, n.sample, burnin,
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=proposal.sd.theta,
                           m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, 
                           target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_loc=target_loc,
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE, self_tune_loc=TRUE,
                           beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial, beta_loc_initial=beta_loc_initial,
                           theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                           prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)
  
  # Check estimated log odds
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  X.standard <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  
  # Save output
  output$description <- paste(sim_name, "_", i, sep="")
  save_output(output, paste("output_", sim_name, "_", i, ".json", sep=""), dst=src)
  
}


###############################
#### Spatial poisson regression
###############################

for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
  
  #########
  ### Cases
  #########
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  d.sub <- d[as.logical(data$locs$status), as.logical(data$locs$status)]
  
  ## Fit model for cases
  # Initial values
  beta_ca_i_ <- rnorm(ncol(X.ca))
  w_i_ <- rnorm(nrow(d.sub))
  phi_i_ <- params$Phi + rnorm(1)
  theta_i_ <- params$Theta + rnorm(1)
  
  # Other parameters
  n.sample_ <- 10000
  burnin_ <- 3000
  L_w_ <- 8
  L_b_ <- 8
  prior_phi_ <- c(3, 40)
  prior_theta_ <- get_gamma_prior(prior_mean=params$Theta, prior_var=5)
  
  # Run fit
  output.sp_ca <- poissonGp(X.ca, Y.ca, d.sub,
                            n.sample=n.sample_, burnin=burnin_, proposal.sd.theta=0.3,
                            L_w=L_w_, L_b=L_b_,
                            beta_initial=beta_ca_i_, w_initial=w_i_, 
                            phi_initial=phi_i_, theta_initial=theta_i_,
                            prior_phi=prior_phi_, prior_theta=prior_theta_)
  
  # Krige random effects
  w.hat_spca <- colMeans(output.sp_ca$samples.w)
  
  kriged_w_ca <- krigeW(output.sp_ca, d, data$locs$ids)
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(data$locs$status))
  
  # Save output and kriged estimates
  output.sp_ca$description <- paste(sim_name, "_", i, sep="")
  save_output(output.sp_ca, paste("output.sp_ca_", sim_name, "_", i, ".json", sep=""), dst=src)
  save_output(kriged_w_ca, paste("output.krige_ca_", sim_name, "_", i, ".json", sep=""), dst=src)
  
  ############
  ### Controls
  ############
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  
  n.sample__ <- 10000
  burnin__ <- 3000
  L_w__ <- 8
  L_b__ <- 8
  
  beta_co_i__ <- rnorm(ncol(X.co))
  w_i__ <- rnorm(nrow(d.sub))
  phi_i__ <- params$Phi + rnorm(1)
  theta_i__ <- params$Theta + rnorm(1)
  
  prior_phi__ <- c(3, 40)
  prior_theta__ <- c(2.5, 2.5)
  
  output.sp_co <- poissonGp(X.co, Y.co, d.sub,
                            n.sample=n.sample__, burnin=burnin__, 
                            L_w=L_w__, L_b=L_b__, proposal.sd.theta=0.3,
                            beta_initial=beta_co_i__, w_initial=w_i__, 
                            phi_initial=phi_i__, theta_initial=theta_i__,
                            prior_phi=prior_phi__, prior_theta=prior_theta__)
  
  # Krige random effects
  w.hat_spco <- colMeans(output.sp_co$samples.w)
  
  kriged_w_co <- krigeW(output.sp_co, d, data$locs$ids)
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(data$locs$status))
  
  output.sp_co$description <- paste(sim_name, "_", i, sep="")
  save_output(output.sp_co, paste("output.sp_co_", sim_name, "_", i, ".json", sep=""), dst=src)
  save_output(kriged_w_co, paste("output.krige_co_", sim_name, "_", i, ".json", sep=""), dst=src)
  
}

####################################
# Bayesian Additive Regression Trees
# (Classifier, High PS Level)
####################################

# Probit BART classifiers are fit to the simulated datasets
# after transformation in the same manner as described in 
# the "low" preferential sampling simulations.

for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  
  #########
  ### Cases
  #########
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  X.coords <- data$locs$coords
  X.combined <- cbind(X.ca, X.coords)
  
  # Generate pseudo-positive points
  Y.new <- c()
  X.new.combined <- c()
  for (idx in 1:length(Y.ca)){
    
    num_ca <- Y.ca[idx]
    Y.new <- c(Y.new, rep(1, num_ca))
    
    # Repeat the covariates for the grid cell for each point
    if (num_ca > 0){
      for (iter in 1:num_ca){
        X.new.combined <- rbind(X.new.combined, X.combined[idx,])
      }
    }
    
  }
  
  ############
  ### Controls
  ############
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  
  for (idx in 1:length(Y.co)){
    
    num_co <- Y.co[idx]
    Y.new <- c(Y.new, rep(0, num_co))
    
    if (num_co > 0){
      for (iter in 1:num_co){
        X.new.combined <- rbind(X.new.combined, X.combined[idx,])
      }
    }
    
  }
  
  #######################
  ### Fit Bart Classifier
  #######################
  
  # By default a probit link function is used if Y has only values 0, 1
  bart_clsf <- bart(X.new.combined, Y.new,
                    # default prior hyperparams from the original BART paper
                    k = 2.0,
                    power = 2.0, 
                    base = 0.95,
                    # as suggested by BART authors
                    ntree = 200,
                    # number of posterior samples
                    ndpost = 2000, 
                    # burnin
                    nskip = 500,
                    printevery = 100, 
                    # no thinning
                    keepevery = 1, 
                    keeptrainfits = TRUE,
                    usequants = FALSE, 
                    numcut = 100, 
                    verbose = TRUE, 
                    # 4 independent MCMC chains
                    nchain = 4, 
                    nthread = 1, 
                    combinechains = TRUE,
                    keeptrees = TRUE)
  
  # Load covariates for all cells in study region
  X.full <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  
  # Combine with cordinates over the entire study region
  r <- caPr.disc[[1]]
  coords_all <- xyFromCell(caPr.disc[[1]], cell = (1:length(r[]))[!is.na(r[])])
  X.full.comb <- cbind(X.full, coords_all)
  
  # Estimate posterior over full study region
  posterior_risk <- predict(bart_clsf, X.full.comb)
  post_mean_risk <- colMeans(posterior_risk)
  
  # Histogram of risk
  hist(post_mean_risk, main = paste("Iteration", i))
  
  ## Save data ------------------------------------------------
  save_output(post_mean_risk, paste("output.bart_clsf_", sim_name, "_", i, ".json", sep=""), dst=src)
  
}
