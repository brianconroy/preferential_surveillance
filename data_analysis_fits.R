
###############################################################
#                       Background                            #
###############################################################

# This script fits a number of models to the plague 
# surveillance dataset (1983-2015), namely:

# (1) The proposed preferential sampling model
#   (1a) Including spatial downscaling
# (2) Spatial poisson models
#   (2a) Including Bayesian kriging
# (3) Probit (BART) classification model

# Risk maps estimated by these models are subsequently 
# compared in the data_analysis_summary.R script, along
# with maps derived from a Poisson Regression model
# (Poisson Regression having been ommitted here due 
# to its simplicity in model fitting). 

library(preferentialSurveillance)

###############################################################
#                       Load Datasets                         #
###############################################################

# PRISM principal components
caPr <- prism_pca
caPr.disc <- aggregate(caPr, fact=5)

# plot Prism PCAS
plot(caPr)

# number of grid cells
N <- n_values(caPr.disc[[1]])
print(N)

# distance matrix of grid cells in study region
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell = cells.all)
d <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))

# mean area of each grid cell
print(mean(raster::area(caPr.disc[[1]])[]))
plot(caPr.disc)

# load processed plague dataset
data <- plague_data

###############################################################
#                   Estimate Initial Values                   #
###############################################################

# Spatial Covariance Function Params --------------------------

# Specify priors for spatial range (theta)
# and marginal variance (phi) of the spatial
# covariance function. These priors hold for
# estimating only the initial spatial random
# effects (below) and not the final model fit.
prior_theta <- get_gamma_prior(prior_mean=3, prior_var=3)
prior_phi <- get_igamma_prior(prior_mean=15, prior_var=3)

# Estimate an initial value for the vector
# of spatial random effects (W). Do so by
# regressing binary indicators for whether
# a grid cell was observed on prism covariates
# plus realizations from a Gaussian process (W)

w_output <- logisticGp(
                       # Binary indicators, covariates, and distance matrix
                       data$loc$status, data$loc$x.scaled, d, 
                       
                       # HMC parameters
                       n.sample=1000, burnin=0, L_beta=8, L_w=8, proposal.sd.theta=0.3,
                       
                       # Initial values
                       w_initial=NULL, theta_initial=NULL, phi_initial=NULL, beta_initial=NULL,
                       
                       # Priors
                       prior_phi=prior_phi, prior_theta=prior_theta
                      )

# optional burnin
w_output <- burnin_logistic_gp_cov(w_output, n.burn=300)

# optionally generate additional MCMC samples
output_ini <- continue_logistic_gp_cov(data, w_output, n.sample=500)

# view traceplots
par(mfrow=c(1,3))
plot(w_output$samples.beta[,1], type='l')
plot(w_output$samples.beta[,2], type='l')
plot(w_output$samples.beta[,3], type='l')

# view acceptance rates
w_output$accept
plot(w_output$samples.w, type='l')
plot(w_output$samples.theta, type='l')
plot(w_output$samples.phi, type='l')

# save the output of initial model fitting
hist(colMeans(w_output$samples.w))
save_output(w_output, "w_inival_output_cdph_v2.json")
# w_output <- load_output("w_inival_output_cdph.json")

# Specify initial values for the locational parameters
# (beta_loc_initial), spatial random effects (w_initial),
# spatial range (theta_initial), and marginal variance
# (phi_initial)
beta_loc_initial <- colMeans(w_output$samples.beta)
w_initial <- colMeans(w_output$samples.w)
theta_initial <- mean(w_output$samples.theta)
phi_initial <- mean(w_output$samples.phi)

# Alpha+ and Beta+ Params -------------------------------------

ini_case <- glm(data$case.data$y ~ data$case.data$x.standardised + w_initial[data$loc$ids] - 1, family='poisson')
alpha_ca_initial <- coefficients(ini_case)[4]
beta_ca_initial <- coefficients(ini_case)[1:3]

# Alpha- and Beta- Params -------------------------------------

ini_ctrl <- glm(data$ctrl.data$y ~ data$ctrl.data$x.standardised + w_initial[data$loc$ids] - 1, family='poisson')
alpha_co_initial <- coefficients(ini_ctrl)[4]
beta_co_initial <- coefficients(ini_ctrl)[1:3]

###############################################################
#                   Fit the Proposed Model                    #
###############################################################

# Specify Prior Distributions ---------------------------------
prior_theta <- get_gamma_prior(prior_mean=theta_initial, prior_var=2)
prior_phi <- get_igamma_prior(prior_mean=phi_initial, prior_var=2)

prior_alpha_ca_mean <- alpha_ca_initial
prior_alpha_co_mean <- alpha_co_initial
prior_alpha_ca_var <- 2
prior_alpha_co_var <- 2

# Specify Hamiltonian Monte Carlo Params Distributions --------

# Number of MCMC samples to generate at first.
# We draw 9,000 more samples subsequently after
# the first 1000. 
n.sample <- 1000

# Initial burnin period. We apply additional burnin
# below (2000 more iterations).
burnin <- 500

# L (HMC parameters): number of steps to simulate
# hamiltonian dynamics within each draw of the samplers.
# Variables are tagged as: 

# _w: spatial random effects
# _ca: covariate params for case counts (beta+)
# _co: covariate params for control counts  (beta-)
# _a_ca: the alpha+ parameter
# _a_co: the alpha- parameter
# _loc: the covariate params for locations (beta loc)

L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
L_loc <- 8

# Proposal standard deviation of the Metropolis-Hastings
# Random Walk to estimate theta (spatial range)
proposal.sd.theta <- 0.15

# Numbers of iterations to self-tune the samplers for 
m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000
m_loc <- 1000

# Target acceptance rates
target_aca=0.65
target_aco=0.65
target_ca=0.65
target_co=0.65
target_w=0.65
target_loc=0.65

# Execute fit -------------------------------------------------
output <- preferentialSampling(data, d, n.sample, burnin, 
                               
                               # HMC params
                               L_w, L_ca, L_co, L_a_ca, L_a_co,
                               proposal.sd.theta=proposal.sd.theta,
                               m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, 
                               target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, 
                               target_co=target_co, target_w=target_w, target_loc=target_loc,
                               self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, 
                               self_tune_ca=TRUE, self_tune_co=TRUE, self_tune_loc=TRUE,
                               
                               # Initial values
                               beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, 
                               alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial, 
                               beta_loc_initial=beta_loc_initial, theta_initial=theta_initial, 
                               phi_initial=phi_initial, w_initial=w_initial,
                               
                               # Prior parametrizations
                               prior_phi=prior_phi, prior_theta=prior_theta, 
                               prior_alpha_ca_var=prior_alpha_ca_var, 
                               prior_alpha_co_var=prior_alpha_co_var)

# Check MCMC acceptance rates
print(output$accept)

# Apply additional burnin, generate more samples --------------

output <- burnin_after(output, n.burn=2000)

output <- continue_mcmc(data, d, output, n.sample=9000)

# View traceplots ---------------------------------------------

# Spatial random effects
plot(apply(output$samples.w, 1, mean), type='l')
view_tr_w(output$samples.w)

# alpha+ and alpha- params
view_tr(output$samples.alpha.ca)
view_tr(output$samples.alpha.co)

mean(output$samples.alpha.ca)
mean(output$samples.alpha.co)

# beta+ params
par(mfrow=c(2,3))
view_tr(output$samples.beta.ca[,1])
view_tr(output$samples.beta.ca[,2])
view_tr(output$samples.beta.ca[,3])

# beta- params
view_tr(output$samples.beta.co[,1])
view_tr(output$samples.beta.co[,2])
view_tr(output$samples.beta.co[,3])

# spatial params
par(mfrow=c(1,1))
view_tr(output$samples.theta)
view_tr(output$samples.phi)

# location related beta params
par(mfrow=c(1,3))
view_tr(output$samples.beta.loc[,1])
view_tr(output$samples.beta.loc[,2])
view_tr(output$samples.beta.loc[,3])

# Save output -------------------------------------------------
output$description <- "v2 cdph baseline "
save_output(output, "output_cdph_baseline_v2.json")

###############################################################
#             Downscale the Spatial Random Effects            #
###############################################################

# We downscale the spatial random effects of the proposed model
# from a coarse to fine grained resolution, following the algorithm
# described in this supplementary webpage:
# http://plague-paper.s3-website-us-east-1.amazonaws.com/

# posterior samples of estimated spatial random effects
samples <- output$samples.w

# define a raster at the desired (16km^2) resolution
r_pred <- caPr[[1]]

# define a raster at the starting (coarse, 426km^2) res.
r_train <- aggregate(r_pred, fact=5)

# interpolate each posterior sample of w via 2d kernel smoothing
w_interp <- interpolate_w_batched(samples, c(0.08, 0.09), r_train, r_pred, batch_size=500)

# save the output: this may take a while as the result is large
save_output(w_interp, "cdph_baseline_interpolated_w_v2.json")

###############################################################
#             Fit the Spatial Poisson Model                   #
###############################################################

## Model for Cases --------------------------------------------

X.ca <- data$case.data$x.standardised
Y.ca <- data$case.data$y
d.sub <- d[as.logical(data$loc$status), as.logical(data$loc$status)]

# Initial parameter values
set.seed(314)
beta_ca_i_ <- rnorm(3)
w_i_ <- rnorm(nrow(d.sub))
phi_i_ <- 10
theta_i_ <- 5

# Hamiltonian Monte Carlo params
n.sample_ <- 10000
burnin_ <- 500
# HMC length params
L_w_ <- 8
L_b_ <- 8

# Prior parametrizations of marginal
# variance and range
prior_phi_ <- c(3, 40)
prior_theta_ <- c(2.5, 2.5)

# Fit the model regressing observed case counts
# on PRISM covariates and spatial random effects
# realized from a Gaussian process
output.sp_ca <- poissonGp(
                          # Covariates, response, distance matrix
                          X.ca, Y.ca, d.sub,
                          
                          # HMC sampling parameters
                          n.sample=n.sample_, burnin=burnin_, proposal.sd.theta=0.3,
                          L_w=L_w_, L_b=L_b_,
                          
                          # Initial values
                          beta_initial=beta_ca_i_, w_initial=w_i_, 
                          phi_initial=phi_i_, theta_initial=theta_i_,
                          
                          # Prior distributions
                          prior_phi=prior_phi_, prior_theta=prior_theta_
                          )

# MCMC acceptance rates
print(output.sp_ca$accept)

# Traceplots
plot(apply(output.sp_ca$samples.w, 1, mean), type='l', col='2')
view_tr_w(output.sp_ca$samples.w)
view_tr(output.sp_ca$samples.beta[,1])
view_tr(output.sp_ca$samples.beta[,2])
view_tr(output.sp_ca$samples.theta)
view_tr(output.sp_ca$samples.phi)
print(colMeans(output.sp_ca$samples.beta))
print(colMeans(output.sp_ca$samples.theta))
print(colMeans(output.sp_ca$samples.phi))

## Krige the estimated spatial random effects -----------------
# Obtain estimates of spatial process at unobserved locations
w.hat_spca <- colMeans(output.sp_ca$samples.w)
beta_ca_sp <- colMeans(output.sp_ca$samples.beta)
kriged_w_ca <- krigeW(output.sp_ca, d, locs$ids)
w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(locs$status))

# Save model output and kriged parameters
output.sp_ca$description <- "cdph_spatial_poisson_case"
save_output(output.sp_ca, "output_cdph_baseline_spatial_poisson_case.json")
save_output(kriged_w_ca, "output_cdph_baseline_krige_ca.json")

## Model for Controls -----------------------------------------

X.co <- data$ctrl.data$x.standardised
Y.co <- data$ctrl.data$y
d.sub <- d[as.logical(data$loc$status), as.logical(data$loc$status)]

# Hamiltonian Monte Carlo params
n.sample__ <- 10000
burnin__ <- 1000
# HMC length params
L_w__ <- 8
L_b__ <- 8

# Initial parameter values
set.seed(314)
beta_co_i__ <- rnorm(3)
w_i__ <- rnorm(nrow(d.sub))
phi_i__ <- 10
theta_i__ <- 5

# Prior parametrizations of marginal
# variance and range
prior_phi__ <- c(3, 40)
prior_theta__ <- c(2.5, 2.5)

# Fit the model regressing observed disease negative counts
# on PRISM covariates and spatial random effects
# realized from a Gaussian process
output.sp_co <- poissonGp(
                          # Covariates, response, distance matrix
                          X.co, Y.co, d.sub,
                          
                          # HMC sampling parameters
                          n.sample=n.sample__, burnin=burnin__, 
                          L_w=L_w__, L_b=L_b__, proposal.sd.theta=0.3,
                          
                          # Initial values
                          beta_initial=beta_co_i__, w_initial=w_i__, 
                          phi_initial=phi_i__, theta_initial=theta_i__,
                          
                          # Prior distributions
                          prior_phi=prior_phi__, prior_theta=prior_theta__
                          )

# MCMC acceptance rates
print(output.sp_co$accept)

# Traceplots
plot(apply(output.sp_co$samples.w, 1, mean), type='l', col='2')
view_tr_w(output.sp_co$samples.w)
view_tr(output.sp_co$samples.beta[,1])
view_tr(output.sp_co$samples.beta[,2])
view_tr(output.sp_co$samples.theta)
view_tr(output.sp_co$samples.phi)
print(colMeans(output.sp_co$samples.beta))
print(colMeans(output.sp_co$samples.theta))
print(colMeans(output.sp_co$samples.phi))

## Krige the estimated spatial random effects -----------------
w.hat_spco <- colMeans(output.sp_co$samples.w)
beta_co_sp <- colMeans(output.sp_co$samples.beta)
kriged_w_co <- krigeW(output.sp_co, d, locs$ids)
w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(locs$status))

# Save model output and kriged parameters
output.sp_co$description <- "cdph_spatial_poisson_ctrl"
save_output(output.sp_co, "output_cdph_baseline_spatial_poisson_ctrl.json")
save_output(kriged_w_co, "output_cdph_baseline_krige_co.json")

###############################################################
#             Fit BART Classification Model                   #
###############################################################

# A probit BART classifier is fit to the (transformed) surveillance
# datasets. We transform our surveillance dataset of aggregated counts 
# into one of binary positive-negative outcomes. Similarly to 
# Carlson et al. (2021), we use the observed Y_i+ positive counts 
# in a grid cell to generate a set of Y_i+ pseudo-positive points 
# located at the centroid of the cell, and Y_i- counts of negative
# specimen to generate a set of Y_i- pseudo-negative points in the 
# same location.

# Here we assume grid cells are of sufficiently small area so that
# the spatial covariates are homogeneous within a given cell, hence
# for simplicity we do not jitter the points to random coordinates
# within the cell. 

# Denoting event of a specimen being positive or negative as B,
# BART estimates Pr(B = 1 | location), the conditional probability 
# a rodent sampled at a location will be plague positive.
# Pr(B = 1 | location) is the quantity of interest, i.e. risk. 

## Cases ------------------------------------------------------------

X.ca <- data$case.data$x.standardised
Y.ca <- data$case.data$y
X.coords <- data$loc$coords
X.combined <- cbind(X.ca, X.coords)
# remove the intercept term
X.combined <- X.combined[,c("x1.standardised", "x2.standardised", "x", "y")]

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

## Controls ---------------------------------------------------------
X.co <- data$ctrl.data$x.standardised
Y.co <- data$ctrl.data$y

# Generate pseudo-negative points
for (idx in 1:length(Y.co)){
  
  num_co <- Y.co[idx]
  Y.new <- c(Y.new, rep(0, num_co))
  
  if (num_co > 0){
    for (iter in 1:num_co){
      X.new.combined <- rbind(X.new.combined, X.combined[idx,])
    }
  }
  
}

## Fit BART Model ---------------------------------------------------------

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

# Load covariates for all cells in study region at fine resolution
X.full <- load_x_ca()
# Rename columns for consistency
colnames(X.full) <- c("", "x1.standardised", "x2.standardised")
# Drop intercept
X.full <- X.full[,c("x1.standardised", "x2.standardised")]

# Combine with cordinates over the entire study region
r <- caPr[[1]]
coords_all <- xyFromCell(r, cell = (1:length(r[]))[!is.na(r[])])
X.full.comb <- cbind(X.full, coords_all)

# Estimate posterior over full study region
posterior_risk <- predict(bart_clsf, X.full.comb)
post_mean_risk <- colMeans(posterior_risk)

# Save output ---------------------------------------------------------
save_output(post_mean_risk, "output_cdph_bart.json")
