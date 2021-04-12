###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (spatial poisson regression,
# poisson regression, and
# BART classifiers) under
# low, high and zero levels of 
# preferential sampling.

# This script evaluates the models
# already fit by simulation.R
###############################

library(preferentialSurveillance)

src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
agg_factor <- 10
n_sims <- 25

#### Prism Principal Components
caPr <- prism_pca
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]


########################################################################
#           Simulation Results: Low Preferential Sampling              #
########################################################################

sampling <- "low"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")

# Define structures to hold metrics
bias_beta_ca_low <- array(NA, c(n_sims, 3))
bias_beta_co_low <- array(NA, c(n_sims, 3))
rmse_ps_low <- c()
rmse_pr_low <- c()
rmse_sp_low <- c()
rmse_bart_low <- c()
mbias_ps_low <- c()
mbias_pr_low <- c()
mbias_sp_low <- c()
mbias_bart_low <- c()
n_cells_obs_low <- c()
prevalences_low <- c()

## Calculate metrics over each simulated dataset -----------------------
for (i in 1:n_sims){
  print(i)
  
  # True parameter values for the iteration
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_", sampling, "_", i, ".json", sep=""), src=src)
  Theta <- params$Theta
  Phi <- params$Phi
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W

  # Preferential sampling model ------------------------------------------
  
  # Load model output
  output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
  
  # Calculate posterior means
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  # Biases
  bias_beta_ca_low[i,] <- beta_ca_h - params$beta.case
  bias_beta_co_low[i,] <- beta_co_h - params$beta.ctrl
  
  # RMSE in log disease odds over all grid cells in the study region
  X.standard <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='Preferential Sampling Model', 
       xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_ps_low <- c(rmse_ps_low, sqrt(mean((lodds.true-lodds.ps)^2)))
  mbias_ps_low <- c(mbias_ps_low, mean(lodds.true-lodds.ps))
  
  # Reference model (Poisson Regression) ---------------------------------
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  # Cases
  rmodel.ca <- glm(Y.ca ~ X.ca-1, family='poisson')
  beta_ca_r <- coefficients(rmodel.ca)
  # Controls
  rmodel.co <- glm(Y.co ~ X.co-1, family='poisson')
  beta_co_r <- coefficients(rmodel.co)
  # Log odds
  lodds.r <- X.standard %*% beta_ca_r - X.standard %*% beta_co_r
  plot(x=lodds.true, y=lodds.r, main='Poisson Regression', 
       xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_pr_low <- c(rmse_pr_low, sqrt(mean((lodds.true-lodds.r)^2)))
  mbias_pr_low <- c(mbias_pr_low, mean(lodds.true-lodds.r))
  
  # Reference model (Spatial Poisson Regression) -------------------------
  output_sp_ca <- load_output(paste("output.sp_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  output_sp_co <- load_output(paste("output.sp_co_", sim_name, "_", i,  ".json", sep=""), src=src)
  kriged_w_ca <- load_output(paste("output.krige_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  kriged_w_co <- load_output(paste("output.krige_co_", sim_name, "_", i, ".json", sep=""), src=src)
  w.hat_spca <- colMeans(output_sp_ca$samples.w)
  w.hat_spco <- colMeans(output_sp_co$samples.w)
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(data$locs$status))
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(data$locs$status))
  beta_co_sp <- colMeans(output_sp_co$samples.beta)
  beta_ca_sp <- colMeans(output_sp_ca$samples.beta)
  lodds_sp <- X.standard %*% beta_ca_sp + w_ca_est - X.standard %*% beta_co_sp - w_co_est
  rmse_sp_low <- c(rmse_sp_low, sqrt(mean((lodds.true-lodds_sp)^2)))
  plot(x=lodds.true, y=lodds_sp, main='Spatial Poisson Regression', 
       xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  mbias_sp_low <- c(mbias_sp_low, mean(lodds.true-lodds_sp))
  
  n_cells_obs_low <- c(n_cells_obs_low, sum(data$locs$status))
  prevalences_low <- c(prevalences_low, sum(data$case.data$y)/(sum(data$case.data$y + data$ctrl.data$y)))
  
  # Probit BART Classification model -------------------------------------
  
  # The BART output is a vector of estimated posterior mean risk values for each
  # grid cell in the study region
  pred_risk_bart <- load_output(paste("output.bart_clsf_", sim_name, "_", i, ".json", sep=""), src=src)

  # Transform posterior mean risk values to obtain the log
  # disease oddds: log[r / (1 - r)] = log odds
  lodds.bart <- log(pred_risk_bart/(1-pred_risk_bart))
  plot(x=lodds.true, y=lodds.bart, main=paste('BART simulation:', i), 
       xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  
  # BART RMSE
  rmse_bart_low <- c(rmse_bart_low, sqrt(mean((lodds.true-lodds.bart)^2)))
  mbias_bart_low <- c(mbias_bart_low, mean(lodds.true-lodds.bart))
  
}

# Initial summary of RMSE in log disease odds across the 4 models
df1 <- data.frame(cbind(rmse_ps_low, rmse_sp_low, rmse_pr_low, rmse_bart_low))
names(df1) <- c("ps", "sp", "pr", "ba")
boxplot(df1, col=rainbow(3, s=0.5))

# Check relationship between RMSE of the proposed model and the 
# number of observed cells in each dataset
plot(x=n_cells_obs_low, y=rmse_ps_low, 
     xlab="Number of Observed Grid Cells", ylab="RMSE (Pref Sampling Model)")


########################################################################
#           Simulation Results: High Preferential Sampling             #
########################################################################


sampling <- "high"
src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")

# Define structures to hold metrics
bias_beta_ca_high <- array(NA, c(n_sims, 3))
bias_beta_co_high <- array(NA, c(n_sims, 3))
rmse_ps_high <- c()
rmse_pr_high <- c()
rmse_sp_high <- c()
rmse_bart_high <- c()
mbias_ps_high <- c()
mbias_pr_high <- c()
mbias_sp_high <- c()
mbias_bart_high <- c()
n_cells_obs_high <- c()
prevalences_high <- c()

## Calculate metrics over each simulated dataset -----------------------
for (i in 1:n_sims){
  print(i)
  
  # True parameter values for the iteration
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
  Theta <- params$Theta
  Phi <- params$Phi
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  # Preferential sampling model ------------------------------------------
  
  # Load model output
  output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
  
  # Calculate posterior means
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  # Biases
  bias_beta_ca_high[i,] <- beta_ca_h - params$beta.case
  bias_beta_co_high[i,] <- beta_co_h - params$beta.ctrl
  
  # RMSE in log disease odds over all grid cells in the study region
  X.standard <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_ps_high <- c(rmse_ps_high, sqrt(mean((lodds.true-lodds.ps)^2)))
  mbias_ps_high <- c(mbias_ps_high, mean(lodds.ps-lodds.true))
  
  # Reference model (Poisson Regression) ---------------------------------
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  # Cases
  rmodel.ca <- glm(Y.ca ~ X.ca-1, family='poisson')
  beta_ca_r <- coefficients(rmodel.ca)
  # Controls
  rmodel.co <- glm(Y.co ~ X.co-1, family='poisson')
  beta_co_r <- coefficients(rmodel.co)
  # Log odds
  lodds.r <- X.standard %*% beta_ca_r - X.standard %*% beta_co_r
  plot(x=lodds.true, y=lodds.r, main='reference', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_pr_high <- c(rmse_pr_high, sqrt(mean((lodds.true-lodds.r)^2)))
  mbias_pr_high <- c(mbias_pr_high, mean(lodds.r-lodds.true))
  
  # Reference model (Spatial Poisson Regression) -------------------------
  output_sp_ca <- load_output(paste("output.sp_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  output_sp_co <- load_output(paste("output.sp_co_", sim_name, "_", i,  ".json", sep=""), src=src)
  kriged_w_ca <- load_output(paste("output.krige_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  kriged_w_co <- load_output(paste("output.krige_co_", sim_name, "_", i, ".json", sep=""), src=src)
  w.hat_spca <- colMeans(output_sp_ca$samples.w)
  w.hat_spco <- colMeans(output_sp_co$samples.w)
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(data$locs$status))
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(data$locs$status))
  beta_co_sp <- colMeans(output_sp_co$samples.beta)
  beta_ca_sp <- colMeans(output_sp_ca$samples.beta)
  lodds_sp <- X.standard %*% beta_ca_sp + w_ca_est - X.standard %*% beta_co_sp - w_co_est
  rmse_sp_high <- c(rmse_sp_high, sqrt(mean((lodds.true-lodds_sp)^2)))
  plot(x=lodds.true, y=lodds_sp, main='reference 2', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  mbias_sp_high <- c(mbias_sp_high, mean(lodds_sp-lodds.true))
  
  n_cells_obs_high <- c(n_cells_obs_high, sum(data$locs$status))
  prevalences_high <- c(prevalences_high, sum(data$case.data$y)/(sum(data$case.data$y + data$ctrl.data$y)))
  
  # Probit BART Classification model -------------------------------------
  
  # The BART output is a vector of estimated posterior mean risk values 
  # for each grid cell in the study region
  pred_risk_bart <- load_output(paste("output.bart_clsf_", sim_name, "_", i, ".json", sep=""), src=src)
  
  # Transform posterior mean risk values to obtain the log disease oddds: 
  # log[r / (1 - r)] = log odds
  lodds.bart <- log(pred_risk_bart/(1-pred_risk_bart))
  plot(x=lodds.true, y=lodds.bart, main=paste('BART simulation:', i), 
       xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  
  # BART RMSE
  rmse_bart_high <- c(rmse_bart_high, sqrt(mean((lodds.true-lodds.bart)^2)))
  mbias_bart_high <- c(mbias_bart_high, mean(lodds.true-lodds.bart))
  
}

# Initial summary of RMSE in log disease odds across the 4 models
boxplot(cbind(rmse_ps_high, rmse_sp_high, rmse_pr_high, rmse_bart_high), col=rainbow(3, s=0.5))

# Check relationship between RMSE of the proposed model and the 
# number of observed cells in each dataset
plot(x=n_cells_obs_high, y=rmse_ps_high, 
     xlab="Number of Observed Grid Cells", ylab="RMSE (Pref Sampling Model)")


########################################################################
#           Simulation Results: No Preferential Sampling             #
########################################################################


sampling <- "none"
src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")

# Define structures to hold metrics
bias_beta_ca_none <- array(NA, c(n_sims, 3))
bias_beta_co_none <- array(NA, c(n_sims, 3))
est_alpha_ca_none <- c()
est_alpha_co_none <- c()
rmse_ps_none <- c()
rmse_pr_none <- c()
rmse_sp_none <- c()
rmse_bart_none <- c()
mbias_ps_none <- c()
mbias_pr_none <- c()
mbias_sp_none <- c()
mbias_bart_none <- c()
n_cells_obs_none <- c()
prevalences_none <- c()

## Calculate metrics over each simulated dataset -----------------------
for (i in 1:n_sims){
  print(i)
  
  # True parameter values for the iteration
  data <- load_output(paste("data_none_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_none_", i, ".json", sep=""), src=src)
  Theta <- params$Theta
  Phi <- params$Phi
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  # Preferential sampling model ------------------------------------------
  
  # Load model output
  output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
  
  # Calculate posterior means
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  # Biases
  bias_beta_ca_none[i,] <- beta_ca_h - params$beta.case
  bias_beta_co_none[i,] <- beta_co_h - params$beta.ctrl
  est_alpha_ca_none <- c(est_alpha_ca_none, alpha_ca_h)
  est_alpha_co_none <- c(est_alpha_co_none, alpha_co_h)
  
  # RMSE in log disease odds over all grid cells in the study region
  X.standard <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_ps_none <- c(rmse_ps_none, sqrt(mean((lodds.true-lodds.ps)^2)))
  mbias_ps_none <- c(mbias_ps_none, mean(lodds.ps-lodds.true))
  
  # Reference model (Poisson Regression) ---------------------------------
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  # Cases
  rmodel.ca <- glm(Y.ca ~ X.ca-1, family='poisson')
  beta_ca_r <- coefficients(rmodel.ca)
  # Controls
  rmodel.co <- glm(Y.co ~ X.co-1, family='poisson')
  beta_co_r <- coefficients(rmodel.co)
  # Log odds
  lodds.r <- X.standard %*% beta_ca_r - X.standard %*% beta_co_r
  plot(x=lodds.true, y=lodds.r, main='reference', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_pr_none <- c(rmse_pr_none, sqrt(mean((lodds.true-lodds.r)^2)))
  mbias_pr_none <- c(mbias_pr_none, mean(lodds.r-lodds.true))
  
  # Reference model (Spatial Poisson Regression) -------------------------
  output_sp_ca <- load_output(paste("output.sp_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  output_sp_co <- load_output(paste("output.sp_co_", sim_name, "_", i,  ".json", sep=""), src=src)
  kriged_w_ca <- load_output(paste("output.krige_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  kriged_w_co <- load_output(paste("output.krige_co_", sim_name, "_", i, ".json", sep=""), src=src)
  w.hat_spca <- colMeans(output_sp_ca$samples.w)
  w.hat_spco <- colMeans(output_sp_co$samples.w)
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(data$locs$status))
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(data$locs$status))
  beta_co_sp <- colMeans(output_sp_co$samples.beta)
  beta_ca_sp <- colMeans(output_sp_ca$samples.beta)
  lodds_sp <- X.standard %*% beta_ca_sp + w_ca_est - X.standard %*% beta_co_sp - w_co_est
  rmse_sp_none <- c(rmse_sp_none, sqrt(mean((lodds.true-lodds_sp)^2)))
  plot(x=lodds.true, y=lodds_sp, main='reference 2', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  mbias_sp_none <- c(mbias_sp_none, mean(lodds_sp-lodds.true))

  n_cells_obs_none <- c(n_cells_obs_none, sum(data$locs$status))
  prevalences_none <- c(prevalences_none, sum(data$case.data$y)/(sum(data$case.data$y + data$ctrl.data$y)))

  # Probit BART Classification model -------------------------------------

  # The BART output is a vector of estimated posterior mean risk values
  # for each grid cell in the study region
  pred_risk_bart <- load_output(paste("output.bart_clsf_", sim_name, "_", i, ".json", sep=""), src=src)

  # Transform posterior mean risk values to obtain the log disease oddds:
  # log[r / (1 - r)] = log odds
  lodds.bart <- log(pred_risk_bart/(1-pred_risk_bart))
  plot(x=lodds.true, y=lodds.bart, main=paste('BART simulation:', i),
       xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

  # BART RMSE
  rmse_bart_none <- c(rmse_bart_none, sqrt(mean((lodds.true-lodds.bart)^2)))
  mbias_bart_none <- c(mbias_bart_none, mean(lodds.true-lodds.bart))
  
}

# Initial summary of RMSE in log disease odds across the 4 models
boxplot(cbind(rmse_ps_none, rmse_sp_none, rmse_pr_none, rmse_bart_none), col=rainbow(3, s=0.5))

# Histograms of estimates for alpha+ and alpha- (should be mean zero)
par(mfrow=c(1,2))
hist(est_alpha_ca_none, main='Alpha(+)', xlab='Posterior Mean'); abline(v=0, col=2)
hist(est_alpha_co_none, main='Alpha(-)', xlab='Posterior Mean'); abline(v=0, col=2)

# Check relationship between RMSE of the proposed model and the 
# number of observed cells in each dataset
plot(x=n_cells_obs_none, y=rmse_ps_none, 
     xlab="Number of Observed Grid Cells", ylab="RMSE (Pref Sampling Model)")


########################################################################
#                     Figure: RMSE Boxplots                            #
########################################################################

## All datasets --------------------------------------------------------

# High preferential sampling
df2 <- data.frame(cbind(c(rmse_ps_high, rmse_sp_high, rmse_pr_high)))#, rmse_bart_high)))
names(df2) <- "RMSE"
df2$Model <- c(rep("Pref-S", length(rmse_ps_high)), 
               rep("Spat-P", length(rmse_ps_high)),
               rep("Pois-R", length(rmse_pr_high)))#,
               #rep("BART", length(rmse_bart_high)))
df2$Sampling <- rep("High", 3*length(rmse_ps_high))

# Low preferential sampling
df1 <- data.frame(cbind(c(rmse_ps_low, rmse_sp_low, rmse_pr_low)))#, rmse_bart_low)))
names(df1) <- "RMSE"
df1$Model <- c(rep("Pref-S", length(rmse_ps_low)), 
               rep("Spat-P", length(rmse_ps_low)),
               rep("Pois-R", length(rmse_ps_low)))#,
               #rep("BART", length(rmse_ps_low)))
df1$Sampling <- rep("Low", 3*length(rmse_ps_low))

# No preferential sampling
df0 <- data.frame(cbind(c(rmse_ps_none, rmse_sp_none, rmse_pr_none)))#, rmse_bart_none)))
names(df0) <- "RMSE"
df0$Model <- c(rep("Pref-S", length(rmse_ps_none)), 
               rep("Spat-P", length(rmse_ps_none)),
               rep("Pois-R", length(rmse_ps_none)))#,
               #rep("BART", length(rmse_ps_none)))
df0$Sampling <- rep("Zero", 3*length(rmse_ps_none))

# Combined
df <- rbind(df0, df1, df2)
df$Model <- factor(df$Model,levels=c("Pref-S", "Spat-P", "Pois-R"))#, "BART"))

# Plot
p1 <- ggplot() + 
  geom_boxplot(data = df, mapping = aes(Sampling, RMSE, fill=Model)) + 
  scale_x_discrete(limits=c("Zero", "Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("A)") + 
  ylim(0,10.5)

## Datasets with >= 75 observed grid cells -----------------------------

high_geq75 <- n_cells_obs_high >= 75
low_geq75 <- n_cells_obs_low >= 75
none_geq75 <- n_cells_obs_none >= 75

# High preferential sampling
df2_75 <- data.frame(cbind(c(rmse_ps_high[high_geq75], 
                             rmse_sp_high[high_geq75], 
                             rmse_pr_high[high_geq75])))#, 
                             #rmse_bart_high[high_geq75])))
names(df2_75) <- "RMSE"
df2_75$Model <- c(rep("Pref-S", sum(high_geq75)), 
                  rep("Spat-P", sum(high_geq75)),
                  rep("Pois-R", sum(high_geq75)))#,
                  #rep("BART", sum(high_geq75)))
df2_75$Sampling <- rep("High", 3*sum(high_geq75))

# Low preferential sampling
df1_75 <- data.frame(cbind(c(rmse_ps_low[low_geq75], 
                             rmse_sp_low[low_geq75], 
                             rmse_pr_low[low_geq75])))#, 
                             #rmse_bart_low[low_geq75])))
names(df1_75) <- "RMSE"
df1_75$Model <- c(rep("Pref-S", sum(low_geq75)), 
                  rep("Spat-P", sum(low_geq75)),
                  rep("Pois-R", sum(low_geq75)))#,
                  #rep("BART", sum(low_geq75)))
df1_75$Sampling <- rep("Low",3*sum(low_geq75))

# No preferential sampling
df0_75 <- data.frame(cbind(c(rmse_ps_none[none_geq75], 
                             rmse_sp_none[none_geq75], 
                             rmse_pr_none[none_geq75])))#, 
                             #rmse_bart_none[none_geq75])))
names(df0_75) <- "RMSE"
df0_75$Model <- c(rep("Pref-S", sum(none_geq75)), 
                  rep("Spat-P", sum(none_geq75)),
                  rep("Pois-R", sum(none_geq75)))#,
                  #rep("BART", sum(none_geq75)))
df0_75$Sampling <- rep("Zero", 3*sum(none_geq75))

# Combined
df_75 <- rbind(df0_75, df1_75, df2_75)
df_75$Model <- factor(df_75$Model,levels=c("Pref-S", "Spat-P", "Pois-R"))#, "BART"))

# Plot
p2 <- ggplot() + 
  geom_boxplot(data = df_75, mapping = aes(Sampling, RMSE, fill=Model)) + 
  scale_x_discrete(limits=c("Zero", "Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("B)") + 
  ylim(0,10.5)

grid.arrange(p1, p2, ncol=2)


########################################################################
#           Figure: RMSE vs observed cells and prevalences             #
########################################################################

# Observed cells
df0 <- data.frame(cbind(n_cells_obs_none, rmse_ps_none))
names(df0) <- c("N", "RMSE")
df0$Sampling <- "Zero"

df1 <- data.frame(cbind(n_cells_obs_high, rmse_ps_high))
names(df1) <- c("N", "RMSE")
df1$Sampling <- "High"

df2 <- data.frame(cbind(n_cells_obs_low, rmse_ps_low))
names(df2) <- c("N", "RMSE")
df2$Sampling <- "Low"

df <- rbind(df0, df1, df2)

p1 <- ggplot(df, aes(x=N, y=RMSE, group=Sampling)) +
  geom_line(aes(color=Sampling)) +
  geom_point(aes(color=Sampling)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C)")

# Prevalences
df_prev1 <- data.frame(cbind(prevalences_high, rmse_ps_high))
names(df_prev1) <- c("Prevalence", "RMSE")
df_prev1$Sampling <- "High"

df_prev2 <- data.frame(cbind(prevalences_low, rmse_ps_low))
names(df_prev2) <- c("Prevalence", "RMSE")
df_prev2$Sampling <- "Low"

df_prev <- rbind(df_prev1, df_prev2)

p2 <- ggplot(df_prev, aes(x=Prevalence, y=RMSE, group=Sampling)) +
  geom_line(aes(color=Sampling)) +
  geom_point(aes(color=Sampling)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("B)")

grid.arrange(p1, p2, ncol=2)


########################################################################
#                 Figure: Study Region Covariates                     #
########################################################################

# Plot the first and second principal components of the PRISM dataset
par(mfrow=c(1,2))
plot(caPr.disc[[1]], main='')
plot(caPr.disc[[2]], main='')


########################################################################
#                       TABLE: RMSE Summary                            #
########################################################################

# Summarize the root mean squared error in estimated log disease oddds
# over all grid cells in the study region. Do so for the 4 models 
# considered in the simulation, under "low" and "high" preferential
# sampling, for all datasets and those with at least 75 observed grid
# cells

rmse_table <- list()
rmse_table[[1]] <- list(
  Sampling="Zero",
  Model="PS",
  Mean_RMSE=round(mean(rmse_ps_none), 3),
  SD_RMSE=round(sd(rmse_ps_none), 3),
  Mean_RMSE_above=round(mean(rmse_ps_none[n_cells_obs_none >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_ps_none[n_cells_obs_none >= 75]), 3)
)
rmse_table[[2]] <- list(
  Sampling="Zero",
  Model="SP",
  Mean_RMSE=round(mean(rmse_sp_none), 3),
  SD_RMSE=round(sd(rmse_sp_none), 3),
  Mean_RMSE_above=round(mean(rmse_sp_none[n_cells_obs_none >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_sp_none[n_cells_obs_none >= 75]), 3)
)
rmse_table[[3]] <- list(
  Sampling="Zero",
  Model="PR",
  Mean_RMSE=round(mean(rmse_pr_none), 3),
  SD_RMSE=round(sd(rmse_pr_none), 3),
  Mean_RMSE_above=round(mean(rmse_pr_none[n_cells_obs_none >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_pr_none[n_cells_obs_none >= 75]), 3)
)
rmse_table[[4]] <- list(
  Sampling="Zero",
  Model="BART",
  Mean_RMSE=round(mean(rmse_bart_none), 3),
  SD_RMSE=round(sd(rmse_bart_none), 3),
  Mean_RMSE_above=round(mean(rmse_bart_none[n_cells_obs_none >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_bart_none[n_cells_obs_none >= 75]), 3)
)
rmse_table[[5]] <- list(
  Sampling="Low",
  Model="PS",
  Mean_RMSE=round(mean(rmse_ps_low), 3),
  SD_RMSE=round(sd(rmse_ps_low), 3),
  Mean_RMSE_above=round(mean(rmse_ps_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_ps_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[6]] <- list(
  Sampling="Low",
  Model="SP",
  Mean_RMSE=round(mean(rmse_sp_low), 3),
  SD_RMSE=round(sd(rmse_sp_low), 3),
  Mean_RMSE_above=round(mean(rmse_sp_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_sp_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[7]] <- list(
  Sampling="Low",
  Model="PR",
  Mean_RMSE=round(mean(rmse_pr_low), 3),
  SD_RMSE=round(sd(rmse_pr_low), 3),
  Mean_RMSE_above=round(mean(rmse_pr_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_pr_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[8]] <- list(
  Sampling="Low",
  Model="BART",
  Mean_RMSE=round(mean(rmse_bart_low), 3),
  SD_RMSE=round(sd(rmse_bart_low), 3),
  Mean_RMSE_above=round(mean(rmse_bart_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_bart_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[9]] <- list(
  Sampling="High",
  Model="PS",
  Mean_RMSE=round(mean(rmse_ps_high), 3),
  SD_RMSE=round(sd(rmse_ps_high), 3),
  Mean_RMSE_above=round(mean(rmse_ps_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_ps_high[n_cells_obs_high >= 75]), 3)
)
rmse_table[[10]] <- list(
  Sampling="High",
  Model="SP",
  Mean_RMSE=round(mean(rmse_sp_high), 3),
  SD_RMSE=round(sd(rmse_sp_high), 3),
  Mean_RMSE_above=round(mean(rmse_sp_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_sp_high[n_cells_obs_high >= 75]), 3)
)
rmse_table[[11]] <- list(
  Sampling="High",
  Model="PR",
  Mean_RMSE=round(mean(rmse_pr_high), 3),
  SD_RMSE=round(sd(rmse_pr_high), 3),
  Mean_RMSE_above=round(mean(rmse_pr_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_pr_high[n_cells_obs_high >= 75]), 3)
)
rmse_table[[12]] <- list(
  Sampling="High",
  Model="BART",
  Mean_RMSE=round(mean(rmse_bart_high), 3),
  SD_RMSE=round(sd(rmse_bart_high), 3),
  Mean_RMSE_above=round(mean(rmse_bart_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_bart_high[n_cells_obs_high >= 75]), 3)
)
write_latex_table(ldply(rmse_table, "data.frame"), fname="rmse_table_bart.txt", 
                  path="/Users/brianconroy/Documents/research/paper_aas/")


########################################################################
#               TABLE: Summary of Simulated Datasets                   #
########################################################################

table_dataset <- list()
table_dataset[[1]] <- list(
  Sampling="Zero",
  Quantity="Observed Cells",
  Mean=round(mean(n_cells_obs_none), 3),
  SD=round(sd(n_cells_obs_none), 3),
  Q1=quantile(n_cells_obs_none, 0.25),
  Q2=quantile(n_cells_obs_none, 0.75)
)
table_dataset[[2]] <- list(
  Sampling="Zero",
  Quantity="Prevalence",
  Mean=round(mean(prevalences_none), 3),
  SD=round(sd(prevalences_none), 3),
  Q1=round(quantile(prevalences_none, 0.25), 3),
  Q2=round(quantile(prevalences_none, 0.75), 3)
)
table_dataset[[3]] <- list(
  Sampling="Low",
  Quantity="Observed Cells",
  Mean=round(mean(n_cells_obs_low), 3),
  SD=round(sd(n_cells_obs_low), 3),
  Q1=quantile(n_cells_obs_low, 0.25),
  Q2=quantile(n_cells_obs_low, 0.75)
)
table_dataset[[4]] <- list(
  Sampling="Low",
  Quantity="Prevalence",
  Mean=round(mean(prevalences_low), 3),
  SD=round(sd(prevalences_low), 3),
  Q1=round(quantile(prevalences_low, 0.25), 3),
  Q2=round(quantile(prevalences_low, 0.75), 3)
)
table_dataset[[5]] <- list(
  Sampling="High",
  Quantity="Observed Cells",
  Mean=round(mean(n_cells_obs_high), 3),
  SD=round(sd(n_cells_obs_high), 3),
  Q1=quantile(n_cells_obs_high, 0.25),
  Q2=quantile(n_cells_obs_high, 0.75)
)
table_dataset[[6]] <- list(
  Sampling="High",
  Quantity="Prevalence",
  Mean=round(mean(prevalences_high), 3),
  SD=round(sd(prevalences_high), 3),
  Q1=round(quantile(prevalences_high, 0.25), 3),
  Q2=round(quantile(prevalences_high, 0.75), 3)
)
write_latex_table(ldply(table_dataset, "data.frame"), fname="dataset_table.txt", 
                  path="/Users/brianconroy/Documents/research/paper_aas/")


########################################################################
#                     ADDITIONAL TABLES & FIGURES                      #
########################################################################

## Biases in Beta parameter estimates ----------------------------------

## Low preferential sampling -------------------------------------------

# Beta case (low preferential sampling)
df_bca_low <- data.frame(c(bias_beta_ca_low[,1], 
                           bias_beta_ca_low[,2], 
                           bias_beta_ca_low[,3]))
names(df_bca_low) <- "Bias"
df_bca_low$Parameter <- c(rep("Beta (+, 0)", n_sims), 
                          rep("Beta (+, 1)", n_sims), 
                          rep("Beta (+, 2)", n_sims))

# Beta ccontrol (low preferential sampling)
df_bco_low <- data.frame(c(bias_beta_co_low[,1], 
                           bias_beta_co_low[,2], 
                           bias_beta_co_low[,3]))
names(df_bco_low) <- "Bias"
df_bco_low$Parameter <- c(rep("Beta (-, 0)", n_sims), 
                          rep("Beta (-, 1)", n_sims), 
                          rep("Beta (-, 2)", n_sims))

df_b_low <- rbind(df_bca_low, df_bco_low)

# Plot of biases
p1 <- ggplot(df_b_low, aes(x=Parameter, y=Bias)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("A)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-12, 10)

## High preferential sampling ------------------------------------------

# Beta case (high preferential sampling)
df_bca_high <- data.frame(c(bias_beta_ca_high[,1], 
                            bias_beta_ca_high[,2], 
                            bias_beta_ca_high[,3]))
names(df_bca_high) <- "Bias"
df_bca_high$Parameter <- c(rep("Beta (+, 0)", n_sims), 
                           rep("Beta (+, 1)", n_sims), 
                           rep("Beta (+, 2)", n_sims))

# Beta control (high preferential sampling)
df_bco_high <- data.frame(c(bias_beta_co_high[,1], 
                            bias_beta_co_high[,2], 
                            bias_beta_co_high[,3]))
names(df_bco_high) <- "Bias"
df_bco_high$Parameter <- c(rep("Beta (-, 0)", n_sims), 
                           rep("Beta (-, 1)", n_sims), 
                           rep("Beta (-, 2)", n_sims))

df_b_high <- rbind(df_bca_high, df_bco_high)

p2 <- ggplot(df_b_high, aes(x=Parameter, y=Bias)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("B)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-12, 10)

## Bias vs number of observed cells ------------------------------------

# Low preferential sampling
df_b_low$N <- c(rep(n_cells_obs_low, 3))

p3 <- ggplot(df_b_low, aes(x=N, y=Bias, group=Parameter)) +
  geom_line(aes(color=Parameter)) +
  geom_point(aes(color=Parameter)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-12, 10)

# High preferential sampling
df_b_high$N <- c(rep(n_cells_obs_high, 3))

p4 <- ggplot(df_b_high, aes(x=N, y=Bias, group=Parameter)) +
  geom_line(aes(color=Parameter)) +
  geom_point(aes(color=Parameter)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("D)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-12, 10)

grid.arrange(p1,p2,p3,p4,ncol=2)

## Boxplots of estimated alpha_+ and alpha_- ---------------------------

# Retrieve parameter estimates -----------------------------------------
# (High preferential sampling)
estimated_aca_high <- c()
estimated_aco_high <- c()
estimated_aca_high_n <- c()
estimated_aco_high_n <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
  output <- load_output(paste("output_", "sim_iteration_v2_high", "_", i, ".json", sep=""), src=src)
  
  aca_hat <- mean(output$samples.alpha.ca)
  aco_hat <- mean(output$samples.alpha.co)
  estimated_aca_high <- c(estimated_aca_high, aca_hat)
  estimated_aco_high <- c(estimated_aco_high, aco_hat)
  
  # if there are over 75 observed cells in the simulated dataset
  if (sum(data$locs$status) >= 75){
    aca_hat <- mean(output$samples.alpha.ca)
    aco_hat <- mean(output$samples.alpha.co)
    estimated_aca_high_n <- c(estimated_aca_high_n, aca_hat)
    estimated_aco_high_n <- c(estimated_aco_high_n, aco_hat)
  }
  
}
estimated_aca_high <- estimated_aca_high - params$alpha.case
estimated_aco_high <- estimated_aco_high - params$alpha.ctrl
estimated_aca_high_n <- estimated_aca_high_n - params$alpha.case
estimated_aco_high_n <- estimated_aco_high_n - params$alpha.ctrl

# Retrieve parameter estimates -----------------------------------------
# (Low preferential sampling)

estimated_aca_low <- c()
estimated_aco_low <- c()
estimated_aca_low_n <- c()
estimated_aco_low_n <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_low_", i, ".json", sep=""), src=src)
  output <- load_output(paste("output_", "sim_iteration_v2_low", "_", i, ".json", sep=""), src=src)
  
  aca_hat <- mean(output$samples.alpha.ca)
  aco_hat <- mean(output$samples.alpha.co)
  estimated_aca_low <- c(estimated_aca_low, aca_hat)
  estimated_aco_low <- c(estimated_aco_low, aco_hat)
  
  # if there are over 75 observed cells in the simulated dataset
  if (sum(data$locs$status) >= 75){
    aca_hat <- mean(output$samples.alpha.ca)
    aco_hat <- mean(output$samples.alpha.co)
    estimated_aca_low_n <- c(estimated_aca_low_n, aca_hat)
    estimated_aco_low_n <- c(estimated_aco_low_n, aco_hat)
  }
  
}
# Calculate biases
estimated_aca_low <- estimated_aca_low - params$alpha.case
estimated_aco_low <- estimated_aco_low - params$alpha.ctrl
estimated_aca_low_n <- estimated_aca_low_n - params$alpha.case
estimated_aco_low_n <- estimated_aco_low_n - params$alpha.ctrl

# Mean biases ------------------------------------------------------------

print(round(mean(estimated_aca_high), 3))
print(round(mean(estimated_aco_high), 3))
print(round(mean(estimated_aca_low), 3))
print(round(mean(estimated_aco_low), 3))

## Create boxplots (All Datasets)  --------------------------------------------------

dfh <- data.frame(c(estimated_aca_high, estimated_aco_high))
names(dfh) <- "Bias"
dfh$Parameter <- c(rep("Alpha (case)", length(estimated_aca_high)), 
                   rep("Alpha (control)", length(estimated_aca_high)))
dfh$Sampling <- rep("High", 2*length(estimated_aca_high))

dfl <- data.frame(c(estimated_aca_low, estimated_aco_low))
names(dfl) <- "Bias"
dfl$Parameter <- c(rep("Alpha (case)", length(estimated_aca_low)), 
                   rep("Alpha (control)", length(estimated_aca_low)))
dfl$Sampling <- rep("Low", 2*length(estimated_aca_low))

dfhl <- rbind(dfh, dfl)

p1 <- ggplot() + 
  geom_boxplot(data = dfhl, mapping = aes(Sampling, Bias, fill=Parameter)) + 
  scale_x_discrete(limits=c("Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("A)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-1, 3)

## Create boxplots (Observed Cells >= 75) --------------------------------------

dfh_n <- data.frame(c(estimated_aca_high_n, estimated_aco_high_n))
names(dfh_n) <- "Bias"
dfh_n$Parameter <- c(rep("Alpha (case)", length(estimated_aca_high_n)), 
                     rep("Alpha (control)", length(estimated_aca_high_n)))
dfh_n$Sampling <- rep("High", 2*length(estimated_aca_high_n))

dfl_n <- data.frame(c(estimated_aca_low_n, estimated_aco_low_n))
names(dfl_n) <- "Bias"
dfl_n$Parameter <- c(rep("Alpha (case)", length(estimated_aca_low_n)), 
                     rep("Alpha (control)", length(estimated_aca_low_n)))
dfl_n$Sampling <- rep("Low", 2*length(estimated_aca_low_n))

dfh_nl <- rbind(dfh_n, dfl_n)

p2 <- ggplot() + 
  geom_boxplot(data = dfh_nl, mapping = aes(Sampling, Bias, fill=Parameter)) + 
  scale_x_discrete(limits=c("Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("B)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-1, 3)
