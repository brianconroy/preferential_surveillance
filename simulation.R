###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (poisson and spatial poisson)
###############################


library(raster)
library(sp)
library(plyr)
library(mvtnorm)
library(MASS)
library(fields)
library(R.utils)
sourceDirectory('R/')


###########################
# Low preferential sampling
###########################
sampling <- "low"
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
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  
  ## Initial Values
  # Random effects, beta.samp, theta, phi
  params <- load_output(paste("params_low_", i, ".json", sep=""), src=src)
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


############################
# High preferential sampling
############################
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
