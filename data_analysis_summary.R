###############################
# Summarizes simulation results
# of the baseline (scurid)
# cdph analysis
###############################

library(plyr)
library(grid)
library(mvtnorm)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('R/')


dst <- "/Users/brianconroy/Documents/research/project1/analysis/"
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=5)


#### Figure: Prism PCAS
plot(caPr)


N <- n_values(caPr.disc[[1]])
print(N)
print(mean(area(caPr.disc[[1]])[]))
plot(caPr.disc)

src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")
output <- load_output("output_cdph_baseline_v2.json")

##################
# data description
##################

coords_all <- cbind(matrix(rodents$Lon_Add_Fix), rodents$Lat_Add_Fix)
loc.disc <- caPr.disc[[1]]
cells_all <- cellFromXY(loc.disc, coords_all)
cells_all <- cells_all[!is.na(cells_all[])]
cells_obs <- sort(unique(cells_all))

loc.disc[][!is.na(loc.disc[])] <- 0
loc.disc[][cells_obs] <- 5

plot(loc.disc)
plot(rasterToPolygons(loc.disc), add=T, border='black', lwd=1) 
points(coords_all, col='2')

##################
# data description
##################

all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
data <- assemble_data(rodents, loc.disc, caPr.disc)
coords <- xyFromCell(caPr.disc, cell=all_ids)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
case.data <- data$case.data
ctrl.data <- data$ctrl.data
locs <- data$loc

# table of summary metrics
count_sums <- rbind(summary(case.data$y), summary(ctrl.data$y))
count_sums <- cbind(c('Plague Positive', 'Plague Negative'), count_sums)
count_sums <- data.frame(count_sums)
count_sums$Total <- c(sum(case.data$y), sum(ctrl.data$y))
names(count_sums) <- c("Disease Status", "Min", "1st Quartile",
                       "Median", "Mean", "3rd Quartile", "Max", "Total")
write_latex_table(count_sums, "cdph_baseline_count_summary.txt", path=dst)

print(sum(case.data$y)/ sum(case.data$y + ctrl.data$y))


###################
# model convergence
###################


# traceplots
fname <- paste("cdph_baseline_traces", ".png", sep="")
png(paste(dst, fname, sep=""),
     width=900, height=700, res=100)
par(mfrow=c(2,3))
plot_traces_general(output)

# parameter estimates
params <- summarize_ps_params(output)
write_latex_table(params, 'cdph_baseline_params.txt', dst)

###########
# risk maps
###########

# random field (downscaled)
w.hat <- colMeans(output$samples.w)
ds <- downscale(w.hat, caPr.disc, caPr)
w.hat_ds <- ds$w.hat_ds
plot(ds$p)

par(mfrow=c(1,2))
rw <- caPr.disc[[1]]
rw[][!is.na(rw[])] <- w.hat
plot(rw)
points(locs$coords, col=1, pch=16)

alpha.ca.hat <- mean(output$samples.alpha.ca)
alpha.co.hat <- mean(output$samples.alpha.co)
beta.ca.hat <- colMeans(output$samples.beta.ca)
beta.co.hat <- colMeans(output$samples.beta.co)
beta.loc.hat <- colMeans(output$samples.beta.loc)

# risk map (low resolution)
X_low <- load_x_ca2(factor=5)
lodds_low <- X_low %*% beta.ca.hat + alpha.ca.hat * w.hat - X_low %*% beta.co.hat - alpha.co.hat * w.hat
risk_low <- calc_risk(lodds_low)

r_lodds_low <- caPr.disc[[1]]
r_lodds_low[][!is.na(r_lodds_low[])] <- lodds_low
plot(r_lodds_low)

r_risk_low <- caPr.disc[[1]]
r_risk_low[][!is.na(r_risk_low[])] <- risk_low
plot(r_risk_low)

# risk map (downscaled)
X_high <- load_x_ca2()
lodds_high <- X_high %*% beta.ca.hat + alpha.ca.hat * w.hat_ds - X_high %*% beta.co.hat - alpha.co.hat * w.hat_ds
risk_high <- calc_risk(lodds_high)

r_lodds_high <- caPr[[2]]
r_lodds_high[][!is.na(r_lodds_high[])] <- lodds_high
plot(r_lodds_high)

r_risk_high <- caPr[[2]]
r_risk_high[][!is.na(r_risk_high[])] <- risk_high
summary(risk_high)


#### Significance map
# perform interpolation
samples <- output$samples.w
r_pred <- caPr[[1]]
r_train <- aggregate(r_pred, fact=5)
w_interp <- interpolate_w_batched(samples, c(0.08, 0.09), r_train, r_pred, batch_size=500)
save_output(w_interp, "cdph_baseline_interpolated_w_v2.json")

w_interp <- load_output("cdph_baseline_interpolated_w_v2.json")

# calculate maps
samples.risk <- calc_posterior_risk_interp(output, w_interp)
sigmaps <- calc_significance(samples.risk, caPr[[1]], threshold=0.064)
# view
plot(sigmaps$r_inds)
plot(sigmaps$r_inds_95)
# sanity check
coords_neg <- cbind(matrix(rodents[rodents$Res=='NEG',]$Lon_Add_Fix), rodents[rodents$Res=='NEG',]$Lat_Add_Fix)
points(coords_neg, col=4, cex=0.2)
coords_pos <- cbind(matrix(rodents[rodents$Res=='POS',]$Lon_Add_Fix), rodents[rodents$Res=='POS',]$Lat_Add_Fix)
points(coords_pos, col=2, cex=0.2)
# export figure
par(mfrow=c(1,2))
plot(sigmaps$r_point_est)
plot(sigmaps$r_inds_95)





#### Figure: covariate contribution to log odds, random field contribution to log odds
X_high <- load_x_ca2()
cov_rodent <- X_high %*% beta.ca.hat - X_high %*% beta.co.hat
w_rodent <- alpha.ca.hat * w.hat_ds - alpha.co.hat * w.hat_ds
r_cov_rodent <- caPr[[1]]
r_cov_rodent[][!is.na(r_cov_rodent[])] <- cov_rodent
r_w_rodent <- caPr[[1]]
r_w_rodent[][!is.na(r_w_rodent[])] <- w_rodent
par(mfrow=c(1,2))
plot(r_cov_rodent, main='A)')
pal <- colorRampPalette(c("blue","red"))
plot(r_w_rodent, main='B)', col=pal(20))
summary(w_rodent[])

#### Figure: risk, contributions
par(mfrow=c(1,2))
plot(r_w_rodent,col=pal(20))
plot(r_cov_rodent)


#################
# Compare to spatial
# poisson model
#################
output.sp_ca <- load_output("output_cdph_baseline_spatial_poisson_case.json")
output.sp_co <- load_output("output_cdph_baseline_spatial_poisson_ctrl.json")
output_krige_ca <- load_output("output_cdph_baseline_krige_ca.json")
output_krige_co <- load_output("output_cdph_baseline_krige_co.json")

w.hat_spca <- colMeans(output.sp_ca$samples.w)
beta_ca_sp <- colMeans(output.sp_ca$samples.beta)
w_ca_est <- combine_w(w.hat_spca, output_krige_ca$mu.new, as.logical(locs$status))

w.hat_spco <- colMeans(output.sp_co$samples.w)
beta_co_sp <- colMeans(output.sp_co$samples.beta)
w_co_est <- combine_w(w.hat_spco, output_krige_co$mu.new, as.logical(locs$status))

X_low <- load_x_ca2(factor=5)
lodds_low_sp <- X_low %*% beta_ca_sp + w_ca_est - X_low %*% beta_co_sp - w_co_est
risk_low_sp <- calc_risk(lodds_low_sp)

plot(x=lodds_low, y=lodds_low_sp); abline(0, 1, col=2)
plot(x=risk_low, y=risk_low_sp); abline(0, 1, col=2)

# downscale
w_ca_ds <- downscale(w_ca_est, caPr.disc, caPr)$w.hat_ds
w_co_ds <- downscale(w_co_est, caPr.disc, caPr)$w.hat_ds

X_high <- load_x_ca2()
lodds_high_sp <- X_high %*% beta_ca_sp + w_ca_ds - X_high %*% beta_co_sp - w_co_ds
risk_high_sp <- calc_risk(lodds_high_sp)

plot(x=lodds_high, y=lodds_high_sp); abline(0, 1, col=2)
plot(x=risk_high, y=risk_high_sp); abline(0, 1, col=2)

r_risk_high_sp <- caPr[[2]]
r_risk_high_sp[][!is.na(r_risk_high_sp[])] <- risk_high_sp


# apply ceiling to estimated values
r_risk_high_sp_c <- r_risk_high_sp
r_risk_high_sp_c[][r_risk_high_sp_c[] > 0.20] = 0.20
plot(r_risk_high_sp_c)
plot(sigmaps$r_point_est)

rescaled <- equalize_scales(sigmaps$r_point_est, r_risk_high_sp_c)
par(mfrow=c(1,2))
r1 <- rescaled[[1]]
r2 <- rescaled[[2]]
r1[][1] <- 0
r2[][1] <- 0
plot(r1)
plot(r2)

plot(rescaled[[1]])
plot(rescaled[[2]])



#### Figure
par(mfrow=c(1,3))
rescaled <- equalize_scales(r_risk_high, r_risk_high_sp)
plot(rescaled[[1]])
plot(rescaled[[2]])
plot(x=risk_high, y=risk_high_sp, xlab='Risk (Proposed Model)', ylab='Risk (Spatial Poisson)'); abline(0, 1, col=2)


#################
# Compare to 
# poisson model
#################

mod.ca <- glm(case.data$y ~ case.data$x.standardised - 1, family='poisson')
mod.co <- glm(ctrl.data$y ~ ctrl.data$x.standardised - 1, family='poisson')

beta.ca.hat_p <- unname(coefficients(mod.ca))
beta.co.hat_p <- unname(coefficients(mod.co))

lodds_low_p <- X_low %*% beta.ca.hat_p - X_low %*% beta.co.hat_p
r_lodds_low_p <- caPr.disc[[1]]
r_lodds_low_p[][!is.na(r_lodds_low_p[])] <- lodds_low_p
plot(r_lodds_low_p)

lodds_high_p <- X_high %*% beta.ca.hat_p - X_high %*% beta.co.hat_p
risk_high_p <- calc_risk(lodds_high_p)

r_lodds_high_p <- caPr[[2]]
r_lodds_high_p[][!is.na(r_lodds_high_p[])] <- lodds_high_p
plot(r_lodds_high_p)

r_risk_high_p <- caPr[[2]]
r_risk_high_p[][!is.na(r_risk_high_p[])] <- risk_high_p

#### Figure
par(mfrow=c(1,2))
rescaled <- equalize_scales(sigmaps$r_point_est, r_risk_high_p)
plot(rescaled[[1]])
plot(rescaled[[2]])
plot(x=risk_high, y=risk_high_p, xlab='Risk (Proposed Model)', ylab='Risk (Poisson)'); abline(0, 1, col=2)


plot(x=risk_high_p, y=risk_high); abline(0,1,col=2)


#########################
# Compare the 3 risk maps
#########################


#### Figure: ps vs poisson risk map
rescaled <- equalize_scales(r_risk_high, r_risk_high_p)
r_risk_high_ <- rescaled[[1]]
r_risk_high_p_ <- rescaled[[2]]

par(mfrow=c(1,2))
plot(r_risk_high_, main='A)')
plot(r_risk_high_p_, main='B)')

diff <- risk_high - risk_high_p
diff <- 5 * (diff > 0)
r_diff <- r_risk_high_p
r_diff[][!is.na(r_diff[])] <- diff
plot(r_diff)


#### Figure: ps vs spatial poisson risk map
rescaled2 <- equalize_scales(r_risk_high, r_risk_high_sp)
r_risk_high_2 <- rescaled2[[1]]
r_risk_high_sp_ <- rescaled2[[2]]

par(mfrow=c(1,2))
plot(r_risk_high_2, main='A)')
plot(r_risk_high_sp_, main='B)')


##################
# Compare per cell 
# risk estimates
##################

par(mfrow=c(1,2))
plot(y=risk_high, 
     x=risk_high_p, main='A)', 
     ylab='Risk (Preferential Sampling)',
     xlab='Risk (Poisson)'); abline(0, 1, col=2)
plot(y=risk_high, 
     x=risk_high_sp, main='B)',
     ylab='Risk (Preferential Sampling)',
     xlab='Risk (Spatial Poisson)'); abline(0, 1, col=2)

# summary table of risk differences
r_diff <- abs(r_risk_high_p[][!is.na(r_risk_high_p[])] - r_risk_high[][!is.na(r_risk_high[])])
r_diff <- 100*r_diff/r_risk_high_p[][!is.na(r_risk_high_p[])]
r_diff_summary <- round(matrix(summary(r_diff), ncol=6), 4)
r_diff_summary <- data.frame(r_diff_summary)

r_diff2 <- abs(r_risk_high_sp[][!is.na(r_risk_high_sp[])] - r_risk_high[][!is.na(r_risk_high[])])
r_diff2 <- 100*r_diff2/r_risk_high_sp[][!is.na(r_risk_high_sp[])]
r_diff_summary2 <- round(matrix(summary(r_diff2), ncol=6), 4)
r_diff_summary2 <- data.frame(r_diff_summary2)

r_diff_summary <- rbind(r_diff_summary, r_diff_summary2)

names(r_diff_summary) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Max")
print(round(r_diff_summary, 3))

# comparison table of parameter estimates
param_comp <- compare_params(output, output.sp_ca, output.sp_co, mod.ca, mod.co)
param_comp$Parameter <- as.character(param_comp$Parameter)
param_comp <- replace_vals(param_comp, 'Parameter', 'Beta 0 (case)', '$\\beta_{0,+}$')
param_comp <- replace_vals(param_comp, 'Parameter', 'Beta 1 (case)', '$\\beta_{1,+}$')
param_comp <- replace_vals(param_comp, 'Parameter', 'Beta 2 (case)', '$\\beta_{2,+}$')
param_comp <- replace_vals(param_comp, 'Parameter', 'Beta 0 (control)', '$\\beta_{0,-}$')
param_comp <- replace_vals(param_comp, 'Parameter', 'Beta 1 (control)', '$\\beta_{1,-}$')
param_comp <- replace_vals(param_comp, 'Parameter', 'Beta 2 (control)', '$\\beta_{2,-}$')
write_latex_table(param_comp, 'cdph_baseline_param_comparison.txt', dst)


#################################
## Inspect regions of high risk
## from the spatial poisson model
#################################


hist(w_ca_est)
hist(w_co_est)

r1 <- caPr.disc[[1]]
r2 <- caPr.disc[[1]]
r1[][!is.na(r1[])] <- w_ca_est
r2[][!is.na(r2[])] <- w_co_est
par(mfrow=c(1,2))
plot(r1, main='A)')
plot(r2, main='B)')

raw_rates <- case.data$y/(case.data$y + ctrl.data$y)
r3 <- caPr.disc[[1]]
r3[][!is.na(r1[])] <- raw_rates
plot(r3)

# compare against rates
w_ca_obs <- w_ca_est[locs$ids]
w_co_obs <- w_co_est[locs$ids]
risk_obs <- risk_low_sp[locs$ids]
y_ca_obs <- case.data$y
risk_rate_wca <- cbind(risk_obs, raw_rates, y_ca_obs, w_ca_obs, w_co_obs)
risk_rate_wca <- data.frame(risk_rate_wca)
plot(x=risk_rate_wca$raw_rates,y=risk_rate_wca$risk_obs); abline(0, 1, col=2)

risk_rate_wca <- risk_rate_wca[with(risk_rate_wca, order(-risk_obs)),]
risk_rate_wca <- round(risk_rate_wca, 3)
head(risk_rate_wca)
write_latex_table(risk_rate_wca[1:5,], 'spatial_poisson_top_cells.txt', dst)

# get corresponding preferential sampling model details
target_ids <- as.numeric(rownames(risk_rate_wca[1:5,]))
raw_rates_target <- raw_rates[target_ids]
y_ca_target <- y_ca_obs[target_ids]
w.hat_target <- w.hat[target_ids]
a_w_ca <- alpha.ca.hat * w.hat_target
a_w_co <- alpha.co.hat * w.hat_target
a_w_difference <- a_w_ca - a_w_co
risk_ps_obs <- risk_low[target_ids]
risk_ps_details <- cbind(risk_ps_obs, raw_rates_target, y_ca_target, a_w_ca, a_w_co, a_w_difference)
risk_ps_details <- data.frame(risk_ps_details)
round(risk_ps_details, 3)
write_latex_table(round(risk_ps_details, 3), 'ps_top_cells_compared.txt', dst)


########################################
## Quantify the strength of preferential 
## sampling apparent in the data
########################################


diffs_ps <- (alpha.ca.hat - alpha.co.hat) * w.hat
diffs_covs <- X_low %*% beta.ca.hat - X_low %*% beta.co.hat
abs_percs <- round(100 * abs(diffs_ps)/(abs(diffs_covs) + abs(diffs_ps)))
hist(abs_percs)
summary(abs_percs)


##########################
## Posterior variance map
##########################

X_rodent <- load_x_standard2(as.logical(data$loc$status), agg_factor=5)
risk_rodent_sep <- calc_posterior_risk(output, X_rodent)
postvar_rodent_sep <- apply(risk_rodent_sep, 2, var)
r <- caPr.disc[[1]]
r[][!is.na(r[])] <- postvar_rodent_sep
plot(r)
