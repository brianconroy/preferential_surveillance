###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (poisson and spatial poisson)
###############################


library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('R/')


sampling <- "none"
prevalence <- "high"
sim_name <- paste("v2_", gen_sim_name(sampling, prevalence), sep="")


#### Simulation parameters
params <- load_output(paste("v2_params_", prevalence, "_", sampling, ".json", sep=""))
Alpha.case <- params$alpha.case
Alpha.ctrl <- params$alpha.ctrl
beta.case <- params$beta.case
beta.ctrl <- params$beta.ctrl
beta.samp <- params$beta.samp
Phi <- params$Phi
Theta <- params$Theta
W <- params$W


#### Simulation data
data <- load_output(paste("v2_data_", prevalence, "_", sampling, ".json", sep=""))
print(sum(data$case.data$y)/sum(data$case.data$y + data$ctrl.data$y))


#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=9)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


################################
#### Preferential sampling model
################################


## Initial Values
# Random effects, beta.samp, theta, phi
prior_theta <- get_gamma_prior(Theta, 5)
prior_phi <- get_igamma_prior(Phi, 5)
w_output <- logisticGpCov(y=data$locs$status, x=data$locs$x.scaled, d, n.sample=1000, burnin=0, L_beta=8, L_w=8, proposal.sd.theta=0.3,
                          w_initial=NULL, theta_initial=NULL, phi_initial=NULL, beta_initial=NULL,
                          prior_phi=prior_phi, prior_theta=prior_theta)

# Optional burnin
w_output <- burnin_logistic_gp_cov(w_output, n.burn=500)

# Optional additional samples
output_ini <- continue_logistic_gp_cov(data, w_output, n.sample=1000)

par(mfrow=c(1,3))
plot(w_output$samples.beta[,1], type='l'); abline(h=beta.samp[1], col=2)
plot(w_output$samples.beta[,2], type='l'); abline(h=beta.samp[2], col=2)
plot(w_output$samples.beta[,3], type='l'); abline(h=beta.samp[3], col=2)
par(mfrow=c(1,1))
plot(x=colMeans(w_output$samples.w), y=W); abline(0, 1, col=2)
plot(w_output$samples.theta, type='l'); abline(h=Theta, col=2)
plot(w_output$samples.phi, type='l'); abline(h=Phi, col=2)

save_output(w_output, "w_inival_output_iteration_v2.json")
# w_output <- load_output("w_inival_output_iteration_v2.json")
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
output <- prefSampleGpV2(data, d, n.sample, burnin,
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=proposal.sd.theta,
                         m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, 
                         target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_loc=target_loc,
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE, self_tune_loc=TRUE,
                         beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial, beta_loc_initial=beta_loc_initial,
                         theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                         prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)

# optional burnin
w_output <- burnin_logistic_gp_cov(w_output, n.burn=300)

# optional additional samples
output_ini <- continue_logistic_gp_cov(data, w_output, n.sample=500)


plot(output$deltas_w)
plot(output$deltas_ca)
plot(output$deltas_co)
plot(output$deltas_aca)
plot(output$deltas_aco)
print(output$accept)

padded_plot(apply(output$samples.w, 1, mean), mean(W))
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output$samples.w, w_true=W)

view_tr(output$samples.theta, Theta)
print(mean(output$samples.theta)); print(Theta)

view_tr(output$samples.phi, Phi)
print(mean(output$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
view_tr(output$samples.beta.ca[,1], beta.case[1])
view_tr(output$samples.beta.ca[,2], beta.case[2])
view_tr(output$samples.beta.ca[,3], beta.case[3])
print(colMeans(output$samples.beta.ca))

padded_plot(output$samples.beta.co[,1], beta.ctrl[1])
padded_plot(output$samples.beta.co[,2], beta.ctrl[2])
padded_plot(output$samples.beta.co[,3], beta.ctrl[3])
print(colMeans(output$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output$samples.alpha.ca, Alpha.case)
view_tr(output$samples.alpha.co, Alpha.ctrl)
print(mean(output$samples.alpha.ca))
print(mean(output$samples.alpha.co))

par(mfrow=c(1,3))
view_tr(output$samples.beta.loc[,1], beta.samp[1])
view_tr(output$samples.beta.loc[,2], beta.samp[2])
view_tr(output$samples.beta.loc[,3], beta.samp[3])

output$description <- sim_name
save_output(output, paste("output_", sim_name, ".json", sep=""))


# w.hat <- colMeans(output$samples.w)
beta_ca_h <- colMeans(output$samples.beta.ca)
beta_co_h <- colMeans(output$samples.beta.co)
alpha_ca_h <- mean(output$samples.alpha.ca)
alpha_co_h <- mean(output$samples.alpha.co)
phi_h <- mean(output$samples.phi)
theta_h <- mean(output$samples.theta)


###############################
#### Spatial poisson regression
###############################


#########
### Cases
#########
X.ca <- data$case.data$x.standardised
Y.ca <- data$case.data$y
d.sub <- d[as.logical(data$locs$status), as.logical(data$locs$status)]

## Fit model for cases
# Initial values
set.seed(314)
beta_ca_i_ <- rnorm(length(beta.case))
w_i_ <- rnorm(nrow(d.sub))
phi_i_ <- Phi + rnorm(1)
theta_i_ <- Theta + rnorm(1)

# Other parameters
n.sample_ <- 10000
burnin_ <- 2000
L_w_ <- 8
L_b_ <- 8
prior_phi_ <- c(3, 40)
prior_theta_ <- get_gamma_prior(prior_mean=Theta, prior_var=5)

# Run fit
output.sp_ca <- poissonGp(X.ca, Y.ca, d.sub,
                          n.sample=n.sample_, burnin=burnin_, proposal.sd.theta=0.3,
                          L_w=L_w_, L_b=L_b_,
                          beta_initial=beta_ca_i_, w_initial=w_i_, 
                          phi_initial=phi_i_, theta_initial=theta_i_,
                          prior_phi=prior_phi_, prior_theta=prior_theta_)

# Optional additional burnin
output.sp_ca <- burnin_poisson_gp(output.sp_ca, n.burn=2000)

# Optional additional samples
output.sp_ca <- continue_poisson_gp(data=data, x=X.ca, y=Y.ca, output_poisson=output.sp_ca, n.sample=200, d.sub=d.sub)

print(output.sp_ca$accept)

# View betas
par(mfrow=c(1,3))
view_tr(output.sp_ca$samples.beta[,1], beta.case[1])
view_tr(output.sp_ca$samples.beta[,2], beta.case[2])
view_tr(output.sp_ca$samples.beta[,3], beta.case[3])
# View random effects
par(mfrow=c(1,1))
plot(apply(output.sp_ca$samples.w, 1, mean), type='l', col='2')
view_tr_w(output.sp_ca$samples.w)
# View covariance parameters
view_tr(output.sp_ca$samples.theta)
view_tr(output.sp_ca$samples.phi)

# Krige random effects
w.hat_spca <- colMeans(output.sp_ca$samples.w)
beta_ca_sp <- colMeans(output.sp_ca$samples.beta)

kriged_w_ca <- krigeW(output.sp_ca, d, locs$ids)
w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(locs$status))

# Save output and kriged estimates
output.sp_ca$description <- paste("spatial_poisson_case", sampling, prevalence, sep="_")
save_output(output.sp_ca, paste("output.sp_ca_", sim_name, ".json", sep=""))
save_params_psc(paste("params.sp_ca_", sim_name, ".json", sep=""))
save_output(kriged_w_ca, paste("output.krige_ca_", sim_name, ".json", sep=""))


############
### Controls
############
X.co <- ctrl.data$x.standardised
Y.co <- ctrl.data$y

n.sample__ <- 16000
burnin__ <- 2000
L_w__ <- 8
L_b__ <- 8

set.seed(314)
beta_co_i__ <- beta.ctrl + rnorm(length(beta.ctrl))
w_i__ <- rnorm(nrow(d.sub))
phi_i__ <- Phi + rnorm(1)
theta_i__ <- Theta + rnorm(1)

prior_phi__ <- c(3, 40)
prior_theta__ <- c(2.5, 2.5)

output.sp_co <- poissonGp(X.co, Y.co, d.sub,
                          n.sample=n.sample__, burnin=burnin__, 
                          L_w=L_w__, L_b=L_b__, proposal.sd.theta=0.3,
                          beta_initial=beta_co_i__, w_initial=w_i__, 
                          phi_initial=phi_i__, theta_initial=theta_i__,
                          prior_phi=prior_phi__, prior_theta=prior_theta__)

print(output.sp_co$accept)

plot(apply(output.sp_co$samples.w, 1, mean), type='l', col='2')
view_tr_w(output.sp_co$samples.w)

view_tr(output.sp_co$samples.beta[,1])
view_tr(output.sp_co$samples.beta[,2])
view_tr(output.sp_co$samples.theta)
view_tr(output.sp_co$samples.phi)
print(colMeans(output.sp_co$samples.beta))
print(colMeans(output.sp_co$samples.theta))
print(colMeans(output.sp_co$samples.phi))

w.hat_spco <- colMeans(output.sp_co$samples.w)
beta_co_sp <- colMeans(output.sp_co$samples.beta)
kriged_w_co <- krigeW(output.sp_co, d, locs$ids)
w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(locs$status))

output.sp_co$description <- paste("spatial_poisson_ctrl", sampling, prevalence, sep="_")
save_output(output.sp_co, paste("output.sp_co_", sim_name, ".json", sep=""))
save_params_psco(paste("params.sp_co_", sim_name, ".json", sep=""))
save_output(kriged_w_co, paste("output.krige_co_", sim_name, ".json", sep=""))


#### Poisson regression
## Cases
rmodel.ca <- glm(Y.ca ~ X.ca-1, family='poisson')
beta_ca_r <- coefficients(rmodel.ca)

## Controls
rmodel.co <- glm(Y.co ~ X.co-1, family='poisson')
beta_co_r <- coefficients(rmodel.co)

save_estimates_pr(beta_ca_r, beta_co_r, paste("estimates_poisson_", sim_name, ".json", sep=""))


####################################
# save and summarize param estimates
# calculate and compare risk surface
####################################


X.standard <- load_x_standard2(as.logical(data$locs$status))
lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
lrisk.true <- lodds.true/(1 - lodds.true)

par(mfrow=c(1, 3))
lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

lodds.sp <- X.standard %*% beta_ca_sp + w_ca_est - X.standard %*% beta_co_sp - w_co_est
lrisk.sp <- lodds.sp/(1-lodds.sp)
plot(x=lodds.true, y=lodds.sp, main='B)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

lodds.r <- X.standard %*% beta_ca_r - X.standard %*% beta_co_r
lrisk.r <- lodds.r/(1-lodds.r)
plot(x=lodds.true, y=lodds.r, main='C)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

surf.true <- overlay(lodds.true, cov.disc)
surf.ps <- overlay(lodds.ps, cov.disc)
surf.sp <- overlay(lodds.sp, cov.disc)
surf.r <- overlay(lodds.r, cov.disc)
par(mfrow=c(2,2))
pal <- colorRampPalette(c("blue","red"))
brk <- seq(-30, 14, by=4)
plot(surf.true, col=pal(12), breaks=brk, main='A)')
plot(surf.ps, col=pal(12), breaks=brk, main='B)')
plot(surf.sp, col=pal(12), breaks=brk, main='C)')
plot(surf.r, col=pal(12), breaks=brk, main='D)')

rsurf.true <- overlay(lrisk.true, cov.disc)
rsurf.ps <- overlay(lrisk.ps, cov.disc)
rsurf.sp <- overlay(lrisk.sp, cov.disc)
rsurf.r <- overlay(lrisk.r, cov.disc)
par(mfrow=c(2,2))
plot(rsurf.true)
plot(rsurf.ps)
plot(rsurf.sp)
plot(rsurf.r)

par(mfrow=c(1,3))
pal <- colorRampPalette(c("green","red"))
diff.ps <- overlay(lodds.ps - lodds.true, cov.disc)
diff.sp <- overlay(lodds.sp - lodds.true, cov.disc)
diff.r <- overlay(lodds.r - lodds.true, cov.disc)
brks <- seq(-25, 15, by=5)
plot(diff.ps, col=pal(9))
plot(diff.sp, col=pal(9))
plot(diff.r, col=pal(9))

rmse.ps <- round(sqrt(mean((lodds.true - lodds.ps)^2)), 2)
rmse.sp <- round(sqrt(mean((lodds.true - lodds.sp)^2)), 2)
rmse.r <- round(sqrt(mean((lodds.true - lodds.r)^2)), 2)
mae.ps <- round(abs(mean(lodds.true - lodds.ps)), 2)
mae.sp <- round(abs(mean(lodds.true - lodds.sp)), 2)
mae.r <- round(abs(mean(lodds.true - lodds.r)), 2)

metrics <- list(
  list(model='preferential sampling', 'rmse'=rmse.ps, 'mae'=mae.ps),
  list(model='spatial regression', 'rmse'=rmse.sp, 'mae'=mae.sp),
  list(model='poisson regression', 'rmse'=rmse.r, 'mae'=mae.r)
)
df <- ldply(metrics, data.frame)
print(df)
