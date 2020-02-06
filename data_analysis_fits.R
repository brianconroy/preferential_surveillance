
library(mvtnorm)
library(R.utils)
sourceDirectory('R/')


caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=5)
N <- n_values(caPr.disc[[1]])
plot(caPr.disc)

src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")

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

# positive counts at each cell
rodents_pos <- rodents[rodents$Res == 'POS',]
coords_pos <- cbind(matrix(rodents_pos$Lon_Add_Fix), rodents_pos$Lat_Add_Fix)
cells_pos <- cellFromXY(loc.disc, coords_pos)
counts_pos <- data.frame(table(cells_pos))
names(counts_pos) <- c('cell', 'count_pos')

# negative counts at each cell
rodents_neg <- rodents[rodents$Res == 'NEG',]
coords_neg <- cbind(matrix(rodents_neg$Lon_Add_Fix), rodents_neg$Lat_Add_Fix)
cells_neg <- cellFromXY(loc.disc, coords_neg)
counts_neg <- data.frame(table(cells_neg))
names(counts_neg) <- c('cell', 'count_neg')

# combine counts
counts_all <- merge(counts_pos, counts_neg, by='cell', all=T)
counts_all$cell <- as.numeric(as.character(counts_all$cell))
counts_all[is.na(counts_all$count_pos),]$count_pos <- 0
counts_all[is.na(counts_all$count_neg),]$count_neg <- 0
counts_all <- counts_all[with(counts_all, order(cell)),]

# location data
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
locs <- list(
  cells=cells_obs,
  status=1 * c(all_ids %in% cells_obs),  
  coords=xyFromCell(loc.disc, cells_obs)
)
locs$ids <- c(1:length(all_ids))[as.logical(locs$status)]
plot(loc.disc)
points(locs$coords)

# covariates
r <- caPr.disc
k <- length(values(r[[1]])[!is.na(values(r[[1]]))])
x <- array(1, c(k, 1))
x.scaled <- x
for (i in 1:length(names(r))){
  x.i <- values(r[[i]])
  x.i <- x.i[!is.na(x.i)]
  x.scaled.i <- (x.i - mean(x.i))/sd(x.i)
  x <- cbind(x, x.i)
  x.scaled <- cbind(x.scaled, x.scaled.i)
}
locs$x.scaled <- x.scaled

# case data
cov.disc <- caPr.disc
x1 <- cov.disc[[1]][][locs$cells]
x2 <- cov.disc[[2]][][locs$cells]
x1.standardised <- (x1 - mean(x1))/sd(x1)
x2.standardised <- (x2 - mean(x2))/sd(x2)
x <- cbind(1, x1, x2)
x.standardised <- cbind(1, x1.standardised, x2.standardised)

case.data <- list(
  y=counts_all$count_pos,
  x.standardised=x.standardised,
  x=x,
  p=3
)
print(sum(case.data$y)) # 1401

# control data
ctrl.data <- list(
  y=counts_all$count_neg,
  x.standardised=x.standardised,
  x=x,
  p=3
)

coords <- xyFromCell(caPr.disc, cell=all_ids)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
data <- list(loc=locs, case.data=case.data, ctrl.data=ctrl.data)

###########
# fit model
###########

prior_theta <- get_gamma_prior(prior_mean=3, prior_var=3)
prior_phi <- get_igamma_prior(prior_mean=15, prior_var=3)

# W initial value
w_output <- logisticGp(locs$status, locs$x.scaled, d, n.sample=1000, burnin=0, L_beta=8, L_w=8, proposal.sd.theta=0.3,
                          w_initial=NULL, theta_initial=NULL, phi_initial=NULL, beta_initial=NULL,
                          prior_phi=prior_phi, prior_theta=prior_theta)

# optional burnin
w_output <- burnin_logistic_gp_cov(w_output, n.burn=300)

# optional additional samples
output_ini <- continue_logistic_gp_cov(data, w_output, n.sample=500)

# view traceplots
par(mfrow=c(1,3))
plot(w_output$samples.beta[,1], type='l')
plot(w_output$samples.beta[,2], type='l')
plot(w_output$samples.beta[,3], type='l')

w_output$accept
view_tr_w(w_output$samples.w)
view_tr(w_output$samples.theta)
view_tr(w_output$samples.phi)

hist(colMeans(w_output$samples.w))
save_output(w_output, "w_inival_output_cdph_v2.json")
# w_output <- load_output("w_inival_output_cdph.json")

beta_loc_initial <- colMeans(w_output$samples.beta)
w_initial <- colMeans(w_output$samples.w)
theta_initial <- mean(w_output$samples.theta)
phi_initial <- mean(w_output$samples.phi)

ini_case <- glm(case.data$y ~ case.data$x.standardised + w_initial[locs$ids] - 1, family='poisson')
alpha_ca_initial <- coefficients(ini_case)[4]
beta_ca_initial <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl.data$y ~ ctrl.data$x.standardised + w_initial[locs$ids] - 1, family='poisson')
alpha_co_initial <- coefficients(ini_ctrl)[4]
beta_co_initial <- coefficients(ini_ctrl)[1:3]

# set prior theta and phi to the values estimated
prior_theta <- get_gamma_prior(prior_mean=theta_initial, prior_var=2)
prior_phi <- get_igamma_prior(prior_mean=phi_initial, prior_var=2)

prior_alpha_ca_mean <- alpha_ca_initial
prior_alpha_co_mean <- alpha_co_initial
prior_alpha_ca_var <- 2
prior_alpha_co_var <- 2

n.sample <- 1000 # 4000
burnin <- 500
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
L_loc <- 8
proposal.sd.theta <- 0.15

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000
m_loc <- 1000
target_aca=0.65
target_aco=0.65
target_ca=0.65
target_co=0.65
target_w=0.65
target_loc=0.65

output <- preferentialSampling(data, d, n.sample, burnin, 
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=proposal.sd.theta,
                         m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, 
                         target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_loc=target_loc,
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE, self_tune_loc=TRUE,
                         beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial, beta_loc_initial=beta_loc_initial,
                         theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                         prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)


print(output$accept)

output <- burnin_after(output, n.burn=2000)

output <- continueMCMC(data, d, output, n.sample=9000)

plot(apply(output$samples.w, 1, mean), type='l')
view_tr_w(output$samples.w)

view_tr(output$samples.alpha.ca)
view_tr(output$samples.alpha.co)

mean(output$samples.alpha.ca)
mean(output$samples.alpha.co)

par(mfrow=c(2,3))
view_tr(output$samples.beta.ca[,1])
view_tr(output$samples.beta.ca[,2])
view_tr(output$samples.beta.ca[,3])

view_tr(output$samples.beta.co[,1])
view_tr(output$samples.beta.co[,2])
view_tr(output$samples.beta.co[,3])

par(mfrow=c(1,1))
view_tr(output$samples.theta)
view_tr(output$samples.phi)

par(mfrow=c(1,3))
view_tr(output$samples.beta.loc[,1])
view_tr(output$samples.beta.loc[,2])
view_tr(output$samples.beta.loc[,3])

output$description <- "v2 cdph baseline "
save_output(output, "output_cdph_baseline_v2.json")


###########################
# Fit spatial poisson model
###########################


X.ca <- case.data$x.standardised
Y.ca <- case.data$y
d.sub <- d[as.logical(locs$status), as.logical(locs$status)]

set.seed(314)
beta_ca_i_ <- rnorm(3)
w_i_ <- rnorm(nrow(d.sub))
phi_i_ <- 10
theta_i_ <- 5

n.sample_ <- 10000
burnin_ <- 500
L_w_ <- 8
L_b_ <- 8

prior_phi_ <- c(3, 40)
prior_theta_ <- c(2.5, 2.5)

output.sp_ca <- poissonGp(X.ca, Y.ca, d.sub,
                          n.sample=n.sample_, burnin=burnin_, proposal.sd.theta=0.3,
                          L_w=L_w_, L_b=L_b_,
                          beta_initial=beta_ca_i_, w_initial=w_i_, 
                          phi_initial=phi_i_, theta_initial=theta_i_,
                          prior_phi=prior_phi_, prior_theta=prior_theta_)

print(output.sp_ca$accept)

plot(apply(output.sp_ca$samples.w, 1, mean), type='l', col='2')
view_tr_w(output.sp_ca$samples.w)

view_tr(output.sp_ca$samples.beta[,1])
view_tr(output.sp_ca$samples.beta[,2])
view_tr(output.sp_ca$samples.theta)
view_tr(output.sp_ca$samples.phi)
print(colMeans(output.sp_ca$samples.beta))
print(colMeans(output.sp_ca$samples.theta))
print(colMeans(output.sp_ca$samples.phi))

w.hat_spca <- colMeans(output.sp_ca$samples.w)
beta_ca_sp <- colMeans(output.sp_ca$samples.beta)
kriged_w_ca <- krigeW(output.sp_ca, d, locs$ids)
w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(locs$status))

output.sp_ca$description <- "cdph_spatial_poisson_case"
save_output(output.sp_ca, "output_cdph_baseline_spatial_poisson_case.json")
save_output(kriged_w_ca, "output_cdph_baseline_krige_ca.json")

## Controls
X.co <- ctrl.data$x.standardised
Y.co <- ctrl.data$y
d.sub <- d[as.logical(locs$status), as.logical(locs$status)]

n.sample__ <- 10000
burnin__ <- 1000
L_w__ <- 8
L_b__ <- 8

set.seed(314)
beta_co_i__ <- rnorm(3)
w_i__ <- rnorm(nrow(d.sub))
phi_i__ <- 10
theta_i__ <- 5

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

output.sp_co$description <- "cdph_spatial_poisson_ctrl"
save_output(output.sp_co, "output_cdph_baseline_spatial_poisson_ctrl.json")
save_output(kriged_w_co, "output_cdph_baseline_krige_co.json")


###########
# downscale
###########

# interpolate w
w.hat <- colMeans(output$samples.w)
rw <- caPr.disc[[1]]
rw[][!is.na(rw[])] <- w.hat

xy <- data.frame(xyFromCell(rw, 1:ncell(rw)))
v <- getValues(rw)

tps <- Tps(xy, v)
p <- raster(caPr[[2]])
p <- interpolate(p, tps)
p <- mask(p, caPr[[1]])
w.hat_ds <- p[][!is.na(p[])]

par(mfrow=c(1,2))
plot(rw)
plot(p)

par(mfrow=c(1,2))
plot(rw)
plot(rw)
points(locs$coords, col=1, pch=16)

################
# calculate risk 
# surfaces
################

alpha.ca.hat <- mean(output$samples.alpha.ca)
alpha.co.hat <- mean(output$samples.alpha.co)
beta.ca.hat <- colMeans(output$samples.beta.ca)
beta.co.hat <- colMeans(output$samples.beta.co)

# low resolution
X_low <- load_x_ca(factor=5)
lodds_low <- X_low %*% beta.ca.hat + alpha.ca.hat * w.hat - X_low %*% beta.co.hat - alpha.co.hat * w.hat
risk_low <- exp(lodds_low/(1 - lodds_low))

r_lodds_low <- caPr.disc[[1]]
r_lodds_low[][!is.na(r_lodds_low[])] <- lodds_low
plot(r_lodds_low)

r_risk_low <- caPr.disc[[1]]
r_risk_low[][!is.na(r_risk_low[])] <- risk_low
plot(r_risk_low)

# high resolution
X_high <- load_x_ca()
lodds_high <- X_high %*% beta.ca.hat + alpha.ca.hat * w.hat_ds - X_high %*% beta.co.hat - alpha.co.hat * w.hat_ds
risk_high <- exp(lodds_high/(1-lodds_high))

r_lodds_high <- caPr[[2]]
r_lodds_high[][!is.na(r_lodds_high[])] <- lodds_high
plot(r_lodds_high)

r_risk_high <- caPr[[2]]
r_risk_high[][!is.na(r_risk_high[])] <- risk_high

pal <- colorRampPalette(c("blue","red"))
plot(r_risk_high, col=pal(15))

par(mfrow=c(1,2))
plot(r_lodds_low)
plot(r_lodds_high)

# compare to naive model
mod.ca <- glm(case.data$y ~ case.data$x.standardised - 1, family='poisson')
mod.co <- glm(ctrl.data$y ~ ctrl.data$x.standardised - 1, family='poisson')

beta.ca.hat_p <- unname(coefficients(mod.ca))
beta.co.hat_p <- unname(coefficients(mod.co))

lodds_low_p <- X_low %*% beta.ca.hat_p - X_low %*% beta.co.hat_p
r_lodds_low_p <- caPr.disc[[1]]
r_lodds_low_p[][!is.na(r_lodds_low_p[])] <- lodds_low_p
plot(r_lodds_low_p)

lodds_high_p <- X_high %*% beta.ca.hat_p - X_high %*% beta.co.hat_p
risk_high_p <- exp(lodds_high_p/(1-lodds_high_p))

r_lodds_high_p <- caPr[[2]]
r_lodds_high_p[][!is.na(r_lodds_high_p[])] <- lodds_high_p
plot(r_lodds_high_p)

r_risk_high_p <- caPr[[2]]
r_risk_high_p[][!is.na(r_risk_high_p[])] <- risk_high_p
plot(r_risk_high_p)

par(mfrow=c(1,2))
plot(r_lodds_high_p)
plot(r_lodds_high)

plot(y=r_lodds_high_p[], x=r_lodds_high[], xlab='log odds (preferential sampling)', ylab='log odds (poisson)'); abline(0, 1, col=2)
plot(y=r_risk_high_p[], x=r_risk_high[], xlab='risk (preferential sampling)', ylab='risk (poisson)'); abline(0, 1, col=2)
