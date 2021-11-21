
###############################################################
#                       Background                            #
###############################################################

# This script compares risk surfaces estimated by the models 
# fit in data_analysis_fits.R. The focus of analysis is the 
# plague surveillance dataset (1983-2015) over the state of 
# California.

# The models compared are:

# (1) The proposed preferential sampling model
#   (1a) Including spatial downscaling
# (2) Spatial poisson models
#   (2a) Including Bayesian kriging
# (3) Probit (BART) classification model
# (4) Poisson regression models

library(preferentialSurveillance)
library(raster)
library(sp)

## Global Variables -------------------------------------------

# Define the destination of summary tables 
# and metrics produced by this script
DST <- "/Users/brianconroy/Documents/research/project1/analysis/"

###############################################################
#                        Load Data                            #
###############################################################

# PRISM principal components
caPr <- prism_pca
caPr.disc <- aggregate(caPr, fact=5)

# plot Prism PCAS
plot(caPr)

# number of grid cells
N <- n_values(caPr.disc[[1]])
print(N)

# mean area of each grid cell
print(mean(raster::area(caPr.disc[[1]])[]))
plot(caPr.disc)

# load processed plague dataset
# Note: this data was removed from the 
# package due to CDPH policy but is 
# available upon request.
rodents <- read.csv(
  "/Users/brianconroy/Documents/research/cdph/data/CDPH_scurid_updated_full.csv",
  header=T,
  sep=",")

# Create the dataset of locations, case and control counts
data <- assemble_data(rodents, caPr.disc)

# extract values
case.data <- data$case.data
ctrl.data <- data$ctrl.data
locs <- data$loc

###############################################################
#                   Load Model Outputs                        #
###############################################################

# load output of the proposed method
output <- load_output("output_cdph_baseline_v2.json")

# Load the spatially downscaled posterior samples of the
# random effects (w). This may take around ~5 min to load.
w_interp <- load_output("cdph_baseline_interpolated_w_v2.json")

# Load outputs of model fitting and kriging
output.sp_ca <- load_output("output_cdph_baseline_spatial_poisson_case.json")
output.sp_co <- load_output("output_cdph_baseline_spatial_poisson_ctrl.json")
output_krige_ca <- load_output("output_cdph_baseline_krige_ca.json")
output_krige_co <- load_output("output_cdph_baseline_krige_co.json")

###############################################################
#                   Study Region Description                  #
###############################################################

# plot the gridded study area with indicators 
# for which cells are observed by the system
plot(rasterToPolygons(data$loc$raster), add=T, border='black', lwd=1)

# plot sampling locations as points
coords_all <- cbind(matrix(rodents$Lon_Add_Fix), rodents$Lat_Add_Fix)
us <- raster::getData("GADM", country="USA", level=1)
ca <- us[match(toupper("California"),toupper(us$NAME_1)),]
plot(ca)
points(coords_all, pch=20, cex=0.5)

## Dataset Description ---------------------------------------

# table of summary metrics
count_sums <- rbind(summary(case.data$y), summary(ctrl.data$y))
count_sums <- cbind(c('Plague Positive', 'Plague Negative'), count_sums)
count_sums <- data.frame(count_sums)
count_sums$Total <- c(sum(case.data$y), sum(ctrl.data$y))
names(count_sums) <- c("Disease Status", "Min", "1st Quartile",
                       "Median", "Mean", "3rd Quartile", "Max", "Total")
write_latex_table(count_sums, "cdph_baseline_count_summary.txt", path=DST)

# Aggregate prevalence of plague over study period
print(sum(case.data$y)/ sum(case.data$y + ctrl.data$y))
# 0.06436349

## Assess Model Convergence ------------------------------------

# Traceplots
fname <- paste("cdph_baseline_traces", ".png", sep="")
png(paste(DST, fname, sep=""),
     width=900, height=700, res=100)
par(mfrow=c(2,3))
plot_traces_general(output)

###############################################################
#                  Proposed Method Risk Map                   #
###############################################################

# Load the covariates (PRISM principal components)
# over the study region at high resolution
X_high <- load_x_ca()

# Calculate posterior risk samples at high resolution
samples.risk <- calc_posterior_risk_ds(output, X_high, w_interp)

# Overlay posterior mean risk onto raster
r_proposed <- caPr[[1]]
r_proposed[][!is.na(r_proposed[])] <- colMeans(samples.risk)
plot(r_proposed)

# Calculate "significance map" showing posterior probabilities
# of risk in each grid cell exceeding 0.064, the statewide
# average risk
sigmaps <- calc_significance(samples.risk, caPr[[1]], threshold=0.064)

# indicators showing cells for which the probability of exceeding
# the statewide average risk is greater than 0.95
plot(sigmaps$r_inds_95)

# Sanity check (overlay cases and controls). This chunk will only
# run if you have acccess to the raw plague surveillance dataset
coords_neg <- cbind(matrix(rodents[rodents$Res=='NEG',]$Lon_Add_Fix), 
                    rodents[rodents$Res=='NEG',]$Lat_Add_Fix)
coords_pos <- cbind(matrix(rodents[rodents$Res=='POS',]$Lon_Add_Fix), 
                    rodents[rodents$Res=='POS',]$Lat_Add_Fix)
points(coords_neg, col=4, cex=0.2)
points(coords_pos, col=2, cex=0.2)

###############################################################
#            Spatial Poisson Regression Risk Map              #
###############################################################

# Parameter estimates ----------------------------------------
w.hat_spca <- colMeans(output.sp_ca$samples.w)
beta_ca_sp <- colMeans(output.sp_ca$samples.beta)
w_ca_est <- combine_w(w.hat_spca, output_krige_ca$mu.new, as.logical(locs$status))

w.hat_spco <- colMeans(output.sp_co$samples.w)
beta_co_sp <- colMeans(output.sp_co$samples.beta)
w_co_est <- combine_w(w.hat_spco, output_krige_co$mu.new, as.logical(locs$status))

# Calculate risk (low resolution) ----------------------------
X_low <- load_x_ca(factor=5)
lodds_low_sp <- X_low %*% beta_ca_sp + w_ca_est - X_low %*% beta_co_sp - w_co_est
risk_low_sp <- calc_risk(lodds_low_sp)

# Downscale --------------------------------------------------

# Downscale the spatial random effects
w_ca_ds <- downscale_tps(w_ca_est, caPr.disc, caPr)$w.hat_ds
w_co_ds <- downscale_tps(w_co_est, caPr.disc, caPr)$w.hat_ds

# High resolution risk maps ----------------------------------
X_high <- load_x_ca()
lodds_high_sp <- X_high %*% beta_ca_sp + w_ca_ds - X_high %*% beta_co_sp - w_co_ds
risk_high_sp <- calc_risk(lodds_high_sp)

r_risk_high_sp <- caPr[[2]]
r_risk_high_sp[][!is.na(r_risk_high_sp[])] <- risk_high_sp

# Apply ceiling to estimated values
r_risk_high_sp_c <- r_risk_high_sp
r_risk_high_sp_c[][r_risk_high_sp_c[] > 0.20] = 0.20
plot(r_risk_high_sp_c)
plot(r_proposed)

# Rescale rasters to have same color scheme --------------------
rescaled <- equalize_scales(r_proposed, r_risk_high_sp_c)
par(mfrow=c(1,2))
r1 <- rescaled[[1]]
r2 <- rescaled[[2]]
r1[][1] <- 0
r2[][1] <- 0
plot(r1)
plot(r2)

###############################################################
#                 Poisson Regression Risk Map                 #
###############################################################

# Fit Poisson regression models -------------------------------

# Fit models regressing case and control counts at **observed**
# grid cells on spatial covariates (PRISM pcas) and intercept
mod.ca <- glm(case.data$y ~ case.data$x.standardised - 1, family='poisson')
mod.co <- glm(ctrl.data$y ~ ctrl.data$x.standardised - 1, family='poisson')

# Estimated coefficients for cases and controls
beta.ca.hat_p <- unname(coefficients(mod.ca))
beta.co.hat_p <- unname(coefficients(mod.co))

# Calculate risk surface (low resolution) ---------------------
lodds_low_p <- X_low %*% beta.ca.hat_p - X_low %*% beta.co.hat_p
r_lodds_low_p <- caPr.disc[[1]]
r_lodds_low_p[][!is.na(r_lodds_low_p[])] <- lodds_low_p
plot(r_lodds_low_p)

# Calculate risk surface (high resolution) --------------------
lodds_high_p <- X_high %*% beta.ca.hat_p - X_high %*% beta.co.hat_p
risk_high_p <- calc_risk(lodds_high_p)

# Plot risk ---------------------------------------------------
r_lodds_high_p <- caPr[[2]]
r_lodds_high_p[][!is.na(r_lodds_high_p[])] <- lodds_high_p
plot(r_lodds_high_p)

r_risk_high_p <- caPr[[2]]
r_risk_high_p[][!is.na(r_risk_high_p[])] <- risk_high_p

# Compare to proposed method ----------------------------------
par(mfrow=c(1,2))
rescaled <- equalize_scales(r_proposed, r_risk_high_p)
plot(rescaled[[1]])
plot(rescaled[[2]])


###############################################################
#                 Model Fit Assessment                        #
###############################################################

# load covariates at coarse resolution
X_low <- load_x_ca(factor=5)

# proposed model -----------------------------------------------------------------------------------------

# parameter estimates from proposed model
beta_ca <- colMeans(output$samples.beta.ca)
beta_co <- colMeans(output$samples.beta.co)
alpha_ca <- mean(output$samples.alpha.ca)
alpha_co <- mean(output$samples.alpha.co)
w <- colMeans(output$samples.w)

# expected disease positive and negative counts at observed grid cells
pos_expected <- exp(X_low %*% beta_ca + alpha_ca * w)[as.logical(data$loc$status)]
neg_expected <- exp(X_low %*% beta_co + alpha_co * w)[as.logical(data$loc$status)]

# poisson model -----------------------------------------------------------------------------------------

# expected disease positive and negative counts at observed grid cells
pos_expected_p <- exp(X_low %*% beta.ca.hat_p)[as.logical(data$loc$status)]
neg_expected_p <- exp(X_low %*% beta.co.hat_p)[as.logical(data$loc$status)]

# spatial poisson model ----------------------------------------------------------------------------------

# expected disease positive and negative counts at observed grid cells
pos_expected_sp <- exp(X_low %*% beta_ca_sp + w_ca_est)[as.logical(data$loc$status)]
neg_expected_sp <- exp(X_low %*% beta_co_sp + w_co_est)[as.logical(data$loc$status)]

# plot and summarize --------------------------------------------------------------------------------------

# observed disease positive and negative counts at observed grid cells
pos_obs <- data$case.data$y
neg_obs <- data$ctrl.data$y

# scatterplots of observed vs expected counts
par(mfrow=c(1,2))
plot(y=pos_expected, x=pos_obs,
     ylab = 'Expected Disease Positive Count',
     xlab = 'Observed Disease Positive Count'
     ); abline(0,1)
plot(y=neg_expected, x=neg_obs,
     ylab = 'Expected Disease Negative Count',
     xlab = 'Observed Disease Negative Count'
); abline(0,1)

plot(y=pos_expected_sp, x=pos_obs,
     ylab = 'Expected Disease Positive Count',
     xlab = 'Observed Disease Positive Count'
); abline(0,1)
plot(y=neg_expected_sp, x=neg_obs,
     ylab = 'Expected Disease Negative Count',
     xlab = 'Observed Disease Negative Count'
); abline(0,1)

# tabularize root mean squared errors 
rmse_fit_tab <- 
  rbind(
    data.frame(
      model = 'proposed',
      rmse_positive = round(sqrt(
        mean((pos_obs - pos_expected)**2)
      ), 3),
      rmse_negative = round(sqrt(
        mean((neg_obs - neg_expected)**2)
      ), 3)
    ), 
    data.frame(
      model = 'spatial poisson',
      rmse_positive = round(sqrt(
        mean((pos_obs - pos_expected_sp)**2)
      ), 3),
      rmse_negative = round(sqrt(
        mean((neg_obs - neg_expected_sp)**2)
      ), 3)
    ),
    data.frame(
      model = 'poisson',
      rmse_positive = round(sqrt(
        mean((pos_obs - pos_expected_p)**2)
      ), 3),
      rmse_negative = round(sqrt(
        mean((neg_obs - neg_expected_p)**2)
      ), 3)
    )
  )

write_latex_table(rmse_fit_tab, "rmse_fit_tab.txt", path=DST)

###############################################################
#                    BART Model Risk Map                      #
###############################################################

# Load posterior mean risk estimates from BART classifier
bart_post_risk <- load_output("output_cdph_bart.json")

# Overlay risk estimates onto raster
r_bart <- caPr[[1]]
r_bart[][!is.na(r_bart[])] <- bart_post_risk

# Plot BART risk map
plot(r_bart)


###############################################################
#                      Additional Plots                       #
###############################################################


# Comparison across thresholds --------------------------------

# To better highlight differences in estimated risk surfaces
# plot areas where the estimated risk of each model exceeds
# a range of thresholds.

plot_indicators <- function(r, threshold, title=''){
  inds <- 1 * (r[][!is.na(r[])] > threshold)
  if (length(unique(inds)) == 1){
    inds[1] <- 0
  }
  r_t <- caPr[[1]]
  r_t[][!is.na(r_t[])] <- inds
  plot(r_t, main=title)
}

# proposed model
par(mfrow=c(2,3))
for (t in c(0.01, 0.02, 0.025, 0.03, 0.033, 0.064)){
  plot_indicators(r_proposed, t, paste('threshold', t))
}
 
# poisson model
par(mfrow=c(3,2))
for (t in c(0.01, 0.02, 0.025, 0.03, 0.033, 0.064)){
  plot_indicators(r_risk_high_p, t, paste('threshold', t))
} 

# spatial poisson model
par(mfrow=c(3,2))
for (t in c(0.01, 0.02, 0.025, 0.03, 0.033, 0.064)){
  plot_indicators(r_risk_high_sp, t, paste('threshold', t))
} 

# Endemic Region Overlay --------------------------------------

# Overlay points of observed cases/controls onto the
# plage "endemic" region, i.e. where the proposed model
# predicts risk to be above a reasonable threshold (taken
# here to be  0.033). 

# NOTE: This section will only run if you
# have access to the raw surveillance data.

# coordinates of cases
rodents_pos <- rodents[rodents$Res == 'POS',]
coords_pos <- cbind(matrix(rodents_pos$Lon_Add_Fix), rodents_pos$Lat_Add_Fix)

# coordinates of controls
rodents_neg <- rodents[rodents$Res == 'NEG',]
coords_neg <- cbind(matrix(rodents_neg$Lon_Add_Fix), rodents_neg$Lat_Add_Fix)

# overlay case and control locations onto "endemic region"
par(mfrow=c(1,1))
plot_indicators(r_proposed, 0.033)
points(coords_neg, col=4, pch=20, cex=0.5)
points(coords_pos, col=2, pch=20, cex=0.5)

###############################################################
#                 Table of Parameter Estimates                #
###############################################################

# load initial values for spatial random effects and covariance params
w_output <- load_output("w_inival_output_cdph_v2.json")

# load initial values for other model parameters
inivals <- load_output("other_initial_values_cdph_v2.json")

summarize_estimate <- function(posterior_samples, param_name, ini_val){
  # summarize the prior and posterior characteristics of an estimate
  return(
    data.frame(
      Parameter = param_name,
      Prior.Mean = round(ini_val, 3),
      Posterior.Mean = round(mean(posterior_samples), 3),
      Posterior.Median = round(median(posterior_samples), 3),
      Credible.Lb = round(quantile(posterior_samples, 0.025), 3),
      Credible.Ub = round(quantile(posterior_samples, 0.975), 3)
    )
  )
}

param_table <- rbind(
  # sampling and spatial params
  summarize_estimate(output$samples.alpha.ca, param_name = '$alpha_+$', ini_val = inivals$alpha_ca_initial),
  summarize_estimate(output$samples.alpha.co, param_name = '$alpha_-$', ini_val = inivals$alpha_co_initial),
  summarize_estimate(output$samples.theta, param_name = '$theta$', ini_val = theta_initial),
  summarize_estimate(output$samples.phi, param_name = '$phi$', ini_val = phi_initial),
  # spatial covariate params (for disease + specimen)
  summarize_estimate(output$samples.beta.ca[,1], param_name = '$beta_+0$', ini_val = inivals$beta_ca_initial[1]),
  summarize_estimate(output$samples.beta.ca[,2], param_name = '$beta_+1$', ini_val = inivals$beta_ca_initial[2]),
  summarize_estimate(output$samples.beta.ca[,3], param_name = '$beta_+2$', ini_val = inivals$beta_ca_initial[3]),
  # spatial covariate params (for disease - specimen)
  summarize_estimate(output$samples.beta.co[,1], param_name = '$beta_-0$', ini_val = inivals$beta_co_initial[1]),
  summarize_estimate(output$samples.beta.co[,2], param_name = '$beta_-1$', ini_val = inivals$beta_co_initial[2]),
  summarize_estimate(output$samples.beta.co[,3], param_name = '$beta_-2$', ini_val = inivals$beta_co_initial[3]),
  # spatial covariate params for the component describing sampling locatioons
  summarize_estimate(output$samples.beta.loc[,1], param_name = '$gamma_0$', ini_val = beta_loc_initial[1]),
  summarize_estimate(output$samples.beta.loc[,2], param_name = '$gamma_1$', ini_val = beta_loc_initial[2]),
  summarize_estimate(output$samples.beta.loc[,3], param_name = '$gamma_2$', ini_val = beta_loc_initial[3])
)

# save in LaTex ready format
write_latex_table(param_table, "param_table.txt", path=DST)

###############################################################
#                 Check Zero Inflation                        #
###############################################################

# percent of disease negative counts equal to zero
print(round(100 * sum(data$ctrl.data$y == 0) / length(data$ctrl.data$y), 3))

# percent of disease positive counts equal to zero
print(round(100 * sum(data$case.data$y == 0) / length(data$case.data$y), 3))
