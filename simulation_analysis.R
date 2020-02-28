###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (poisson and spatial poisson)
###############################


library(raster)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('R/')


src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
agg_factor <- 10
n_sims <- 25

#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


######
# RMSE: low 
######


sampling <- "low"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")

bias_beta_ca_low <- array(NA, c(n_sims, 3))
bias_beta_co_low <- array(NA, c(n_sims, 3))
rmse_ps_low <- c()
rmse_pr_low <- c()
rmse_sp_low <- c()
mbias_ps_low <- c()
mbias_pr_low <- c()
mbias_sp_low <- c()
n_cells_obs_low <- c()
prevalences_low <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_", sampling, "_", i, ".json", sep=""), src=src)
  Theta <- params$Theta
  Phi <- params$Phi
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
  
  # Check estimated log odds
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  bias_beta_ca_low[i,] <- beta_ca_h - params$beta.case
  bias_beta_co_low[i,] <- beta_co_h - params$beta.ctrl
  
  # RMSE
  X.standard <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_ps_low <- c(rmse_ps_low, sqrt(mean((lodds.true-lodds.ps)^2)))
  mbias_ps_low <- c(mbias_ps_low, mean(lodds.true-lodds.ps))
  
  # Reference model (Poisson Regression)
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
  rmse_pr_low <- c(rmse_pr_low, sqrt(mean((lodds.true-lodds.r)^2)))
  mbias_pr_low <- c(mbias_pr_low, mean(lodds.true-lodds.r))
  
  # Reference model (Spatial Poisson Regression)
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
  plot(x=lodds.true, y=lodds_sp, main='reference 2', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  mbias_sp_low <- c(mbias_sp_low, mean(lodds.true-lodds_sp))
  
  n_cells_obs_low <- c(n_cells_obs_low, sum(data$locs$status))
  prevalences_low <- c(prevalences_low, sum(data$case.data$y)/(sum(data$case.data$y + data$ctrl.data$y)))
  
}


summary(rmse_ps_low)
summary(rmse_sp_low)
summary(rmse_pr_low)

df1 <- data.frame(cbind(rmse_ps_low, rmse_sp_low, rmse_pr_low))
names(df1) <- c("ps", "sp", "pr")
boxplot(df1, col=rainbow(3, s=0.5))
plot(x=n_cells_obs_low, y=rmse_ps_low)
points(x=n_cells_obs_low, y=rmse_pr_low, col=2)


######
# RMSE: high
######


sampling <- "high"
src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")

bias_beta_ca_high <- array(NA, c(n_sims, 3))
bias_beta_co_high <- array(NA, c(n_sims, 3))
rmse_ps_high <- c()
rmse_pr_high <- c()
rmse_sp_high <- c()
mbias_ps_high <- c()
mbias_pr_high <- c()
mbias_sp_high <- c()
n_cells_obs_high <- c()
prevalences_high <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
  Theta <- params$Theta
  Phi <- params$Phi
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
  
  # Check estimated log odds
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  bias_beta_ca_high[i,] <- beta_ca_h - params$beta.case
  bias_beta_co_high[i,] <- beta_co_h - params$beta.ctrl
  
  # # RMSE
  X.standard <- load_x_standard(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_ps_high <- c(rmse_ps_high, sqrt(mean((lodds.true-lodds.ps)^2)))
  mbias_ps_high <- c(mbias_ps_high, mean(lodds.ps-lodds.true))
  
  # Reference model (Poisson Regression)
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
  
  # Reference model (Spatial Poisson Regression)
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
  
}


# initial inspection
summary(rmse_ps_high)
summary(rmse_sp_high)
summary(rmse_pr_high)

boxplot(cbind(rmse_ps_high, rmse_sp_high, rmse_pr_high), col=rainbow(3, s=0.5))
plot(x=n_cells_obs_high, y=rmse_ps_high)
points(x=n_cells_obs_high, y=rmse_pr_high, col=2)


#########
# Figure: RMSE boxplots
#########

## All datasets

# High
df2 <- data.frame(cbind(c(rmse_ps_high, rmse_sp_high, rmse_pr_high)))
names(df2) <- "RMSE"
df2$Model <- c(rep("PS", length(rmse_ps_high)), 
               rep("SP", length(rmse_ps_high)),
               rep("PR", length(rmse_pr_high)))
df2$Sampling <- rep("High", 3*length(rmse_ps_high))

# Low
df1 <- data.frame(cbind(c(rmse_ps_low, rmse_sp_low, rmse_pr_low)))
names(df1) <- "RMSE"
df1$Model <- c(rep("PS", length(rmse_ps_high)), 
               rep("SP", length(rmse_ps_high)),
               rep("PR", length(rmse_pr_high)))
df1$Sampling <- rep("Low", 3*length(rmse_ps_high))

# Combined
df <- rbind(df1, df2)
df$Model <- factor(df$Model,levels=c("PS", "SP", "PR"))

# Plot
p1 <- ggplot() + 
  geom_boxplot(data = df, mapping = aes(Sampling, RMSE, fill=Model)) + 
  scale_x_discrete(limits=c("Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("A)") + 
  ylim(0,10.5)

## Datasets with >= 75 observed cells

high_geq75 <- n_cells_obs_high >= 75
low_geq75 <- n_cells_obs_low >= 75

# High
df2_75 <- data.frame(cbind(c(rmse_ps_high[high_geq75], rmse_sp_high[high_geq75], rmse_pr_high[high_geq75])))
names(df2_75) <- "RMSE"
df2_75$Model <- c(rep("PS", sum(high_geq75)), 
                  rep("SP", sum(high_geq75)),
                  rep("PR", sum(high_geq75)))
df2_75$Sampling <- rep("High", 3*sum(high_geq75))

# Low
df1_75 <- data.frame(cbind(c(rmse_ps_low[low_geq75], rmse_sp_low[low_geq75], rmse_pr_low[low_geq75])))
names(df1_75) <- "RMSE"
df1_75$Model <- c(rep("PS", sum(low_geq75)), 
                  rep("SP", sum(low_geq75)),
                  rep("PR", sum(low_geq75)))
df1_75$Sampling <- rep("Low", 3*sum(low_geq75))

# Combined
df_75 <- rbind(df1_75, df2_75)
df_75$Model <- factor(df_75$Model,levels=c("PS", "SP", "PR"))

# Plot
p2 <- ggplot() + 
  geom_boxplot(data = df_75, mapping = aes(Sampling, RMSE, fill=Model)) + 
  scale_x_discrete(limits=c("Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("B)") + 
  ylim(0,10.5)

grid.arrange(p1, p2, ncol=2)

#########
# Figure: RMSE vs observed cells and prevalences
#########

# Observed cells
df1 <- data.frame(cbind(n_cells_obs_high, rmse_ps_high))
names(df1) <- c("N", "RMSE")
df1$Sampling <- "High"

df2 <- data.frame(cbind(n_cells_obs_low, rmse_ps_low))
names(df2) <- c("N", "RMSE")
df2$Sampling <- "Low"

df <- rbind(df1, df2)

p1 <- ggplot(df, aes(x=N, y=RMSE, group=Sampling)) +
  geom_line(aes(color=Sampling)) +
  geom_point(aes(color=Sampling)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("A)")

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


#########
# Figure: Plot of study region covariates
#########
par(mfrow=c(1,2))
plot(caPr.disc[[1]], main='')
plot(caPr.disc[[2]], main='')


########
# Table: RMSE summary
########
rmse_table <- list()
rmse_table[[1]] <- list(
  Sampling="Low",
  Model="PS",
  Mean_RMSE=round(mean(rmse_ps_low), 3),
  SD_RMSE=round(sd(rmse_ps_low), 3),
  Mean_RMSE_above=round(mean(rmse_ps_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_ps_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[2]] <- list(
  Sampling="Low",
  Model="SP",
  Mean_RMSE=round(mean(rmse_sp_low), 3),
  SD_RMSE=round(sd(rmse_sp_low), 3),
  Mean_RMSE_above=round(mean(rmse_sp_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_sp_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[3]] <- list(
  Sampling="Low",
  Model="PR",
  Mean_RMSE=round(mean(rmse_pr_low), 3),
  SD_RMSE=round(sd(rmse_pr_low), 3),
  Mean_RMSE_above=round(mean(rmse_pr_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_pr_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[4]] <- list(
  Sampling="High",
  Model="PS",
  Mean_RMSE=round(mean(rmse_ps_high), 3),
  SD_RMSE=round(sd(rmse_ps_high), 3),
  Mean_RMSE_above=round(mean(rmse_ps_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_ps_high[n_cells_obs_high >= 75]), 3)
)
rmse_table[[5]] <- list(
  Sampling="High",
  Model="SP",
  Mean_RMSE=round(mean(rmse_sp_high), 3),
  SD_RMSE=round(sd(rmse_sp_high), 3),
  Mean_RMSE_above=round(mean(rmse_sp_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_sp_high[n_cells_obs_high >= 75]), 3)
)
rmse_table[[6]] <- list(
  Sampling="High",
  Model="PR",
  Mean_RMSE=round(mean(rmse_pr_high), 3),
  SD_RMSE=round(sd(rmse_pr_high), 3),
  Mean_RMSE_above=round(mean(rmse_pr_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_pr_high[n_cells_obs_high >= 75]), 3)
)
write_latex_table(ldply(rmse_table, "data.frame"), fname="rmse_table.txt", 
                  path="/Users/brianconroy/Documents/research/paper/")


########
# Table: Dataset summary
########
table_dataset <- list()
table_dataset[[1]] <- list(
  Sampling="Low",
  Quantity="Observed Cells",
  Mean=round(mean(n_cells_obs_low), 3),
  SD=round(sd(n_cells_obs_low), 3),
  Q1=quantile(n_cells_obs_low, 0.25),
  Q2=quantile(n_cells_obs_low, 0.75)
)
table_dataset[[2]] <- list(
  Sampling="Low",
  Quantity="Prevalence",
  Mean=round(mean(prevalences_low), 3),
  SD=round(sd(prevalences_low), 3),
  Q1=round(quantile(prevalences_low, 0.25), 3),
  Q2=round(quantile(prevalences_low, 0.75), 3)
)
table_dataset[[3]] <- list(
  Sampling="High",
  Quantity="Observed Cells",
  Mean=round(mean(n_cells_obs_high), 3),
  SD=round(sd(n_cells_obs_high), 3),
  Q1=quantile(n_cells_obs_high, 0.25),
  Q2=quantile(n_cells_obs_high, 0.75)
)
table_dataset[[4]] <- list(
  Sampling="High",
  Quantity="Prevalence",
  Mean=round(mean(prevalences_high), 3),
  SD=round(sd(prevalences_high), 3),
  Q1=round(quantile(prevalences_high, 0.25), 3),
  Q2=round(quantile(prevalences_high, 0.75), 3)
)
write_latex_table(ldply(table_dataset, "data.frame"), fname="dataset_table.txt", 
                  path="/Users/brianconroy/Documents/research/paper/")

# Biases in beta
# plot 1, low. plot 2, high
df_bca_low <- data.frame(c(bias_beta_ca_low[,1], bias_beta_ca_low[,2], bias_beta_ca_low[,3]))
names(df_bca_low) <- "Bias"
df_bca_low$Parameter <- c(rep("Beta (+, 0)", n_sims), rep("Beta (+, 1)", n_sims), rep("Beta (+, 2)", n_sims))

df_bco_low <- data.frame(c(bias_beta_co_low[,1], bias_beta_co_low[,2], bias_beta_co_low[,3]))
names(df_bco_low) <- "Bias"
df_bco_low$Parameter <- c(rep("Beta (-, 0)", n_sims), rep("Beta (-, 1)", n_sims), rep("Beta (-, 2)", n_sims))

df_b_low <- rbind(df_bca_low, df_bco_low)

p1 <- ggplot(df_b_low, aes(x=Parameter, y=Bias)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("A)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-12, 10)

# High
df_bca_high <- data.frame(c(bias_beta_ca_high[,1], bias_beta_ca_high[,2], bias_beta_ca_high[,3]))
names(df_bca_high) <- "Bias"
df_bca_high$Parameter <- c(rep("Beta (+, 0)", n_sims), rep("Beta (+, 1)", n_sims), rep("Beta (+, 2)", n_sims))

df_bco_high <- data.frame(c(bias_beta_co_high[,1], bias_beta_co_high[,2], bias_beta_co_high[,3]))
names(df_bco_high) <- "Bias"
df_bco_high$Parameter <- c(rep("Beta (-, 0)", n_sims), rep("Beta (-, 1)", n_sims), rep("Beta (-, 2)", n_sims))

df_b_high <- rbind(df_bca_high, df_bco_high)

p2 <- ggplot(df_b_high, aes(x=Parameter, y=Bias)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("B)") + 
  geom_hline(yintercept=0, color="red") + 
  ylim(-12, 10)

# Bias vs sample size
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

# High
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

########
# Figure: Boxplots of estimated alpha_+ and alpha_-
########

## All datasets

# High prevalence
estimated_aca_high <- c()
estimated_aco_high <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
  output <- load_output(paste("output_", "sim_iteration_v2_high", "_", i, ".json", sep=""), src=src)
  
  aca_hat <- mean(output$samples.alpha.ca)
  aco_hat <- mean(output$samples.alpha.co)
  estimated_aca_high <- c(estimated_aca_high, aca_hat)
  estimated_aco_high <- c(estimated_aco_high, aco_hat)
  
}
estimated_aca_high <- estimated_aca_high - params$alpha.case
estimated_aco_high <- estimated_aco_high - params$alpha.ctrl

# Low prevalence
estimated_aca_low <- c()
estimated_aco_low <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_low_", i, ".json", sep=""), src=src)
  output <- load_output(paste("output_", "sim_iteration_v2_low", "_", i, ".json", sep=""), src=src)
  
  aca_hat <- mean(output$samples.alpha.ca)
  aco_hat <- mean(output$samples.alpha.co)
  estimated_aca_low <- c(estimated_aca_low, aca_hat)
  estimated_aco_low <- c(estimated_aco_low, aco_hat)
  
}
estimated_aca_low <- estimated_aca_low - params$alpha.case
estimated_aco_low <- estimated_aco_low - params$alpha.ctrl

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

print(round(mean(estimated_aca_high), 3))
print(round(mean(estimated_aco_high), 3))
print(round(mean(estimated_aca_low), 3))
print(round(mean(estimated_aco_low), 3))

## n >= 75 datasets

# High prevalence
estimated_aca_high_n <- c()
estimated_aco_high_n <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
  
  if (sum(data$locs$status) >= 75){
    output <- load_output(paste("output_", "sim_iteration_v2_high", "_", i, ".json", sep=""), src=src)
    
    aca_hat <- mean(output$samples.alpha.ca)
    aco_hat <- mean(output$samples.alpha.co)
    estimated_aca_high_n <- c(estimated_aca_high_n, aca_hat)
    estimated_aco_high_n <- c(estimated_aco_high_n, aco_hat)
  }
  
}
estimated_aca_high_n <- estimated_aca_high_n - params$alpha.case
estimated_aco_high_n <- estimated_aco_high_n - params$alpha.ctrl

estimated_aca_low_n <- c()
estimated_aco_low_n <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_low_", i, ".json", sep=""), src=src)
  
  if (sum(data$locs$status) >= 75){
    output <- load_output(paste("output_", "sim_iteration_v2_low", "_", i, ".json", sep=""), src=src)
    
    aca_hat <- mean(output$samples.alpha.ca)
    aco_hat <- mean(output$samples.alpha.co)
    estimated_aca_low_n <- c(estimated_aca_low_n, aca_hat)
    estimated_aco_low_n <- c(estimated_aco_low_n, aco_hat)
  }
  
}
estimated_aca_low_n <- estimated_aca_low_n - params$alpha.case
estimated_aco_low_n <- estimated_aco_low_n - params$alpha.ctrl

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

# Plot biases vs sample size
dfnb_low <- data.frame(cbind(n_cells_obs_low, estimated_aca_low))
names(dfnb_low) <- c("N", "Bias")
dfnb_low$Parameter <- "Alpha (case)"

dfnb_low_co <- data.frame(cbind(n_cells_obs_low, estimated_aco_low))
names(dfnb_low_co) <- c("N", "Bias")
dfnb_low_co$Parameter <- "Alpha (control)"

dfnb_low <- rbind(dfnb_low, dfnb_low_co)

p3 <- ggplot(dfnb_low, aes(x=N, y=Bias, group=Parameter)) +
  geom_line(aes(color=Parameter)) +
  geom_point(aes(color=Parameter)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C)") + 
  ylim(-1, 3) + 
  geom_hline(yintercept=0, color="red")


n_cells_obs_low <- c()
for (i in 1:n_sims){
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  n_cells_obs_low <- c(n_cells_obs_low, sum(data$locs$status))
}

n_cells_obs_high <- c()
for (i in 1:n_sims){
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  n_cells_obs_high <- c(n_cells_obs_high, sum(data$locs$status))
}

dfnb_high <- data.frame(cbind(n_cells_obs_high, estimated_aca_high))
names(dfnb_high) <- c("N", "Bias")
dfnb_high$Parameter <- "Alpha (case)"

dfnb_high_co <- data.frame(cbind(n_cells_obs_high, estimated_aco_high))
names(dfnb_high_co) <- c("N", "Bias")
dfnb_high_co$Parameter <- "Alpha (control)"

dfnb_high <- rbind(dfnb_high, dfnb_high_co)

p4 <- ggplot(dfnb_high, aes(x=N, y=Bias, group=Parameter)) +
  geom_line(aes(color=Parameter)) +
  geom_point(aes(color=Parameter)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("D)") + 
  ylim(-1, 3) + 
  geom_hline(yintercept=0, color="red")

grid.arrange(p1, p2, p3, p4, ncol=2)
