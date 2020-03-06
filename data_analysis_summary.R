###############################
# Summarizes simulation results
# of the baseline (scurid)
# cdph analysis
###############################

library(raster)
library(plyr)
library(grid)
library(mvtnorm)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('R/')


dst <- "/Users/brianconroy/Documents/research/project1/analysis/"
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=5)


#### Figure: Prism PCAS
plot(caPr)

N <- n_values(caPr.disc[[1]])
print(N)
print(mean(raster::area(caPr.disc[[1]])[]))
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

# perform downscaling
samples <- output$samples.w
r_pred <- caPr[[1]]
r_train <- aggregate(r_pred, fact=5)
# w_interp <- interpolate_w_batched(samples, c(0.08, 0.09), r_train, r_pred, batch_size=500)
# save_output(w_interp, "cdph_baseline_interpolated_w_v2.json")
w_interp <- load_output("cdph_baseline_interpolated_w_v2.json")

# calculate maps
samples.risk <- calc_posterior_risk(output, w_interp)
sigmaps <- calc_significance(samples.risk, caPr[[1]], threshold=0.064)

# view
plot(sigmaps$r_inds)
plot(sigmaps$r_inds_95)

# sanity check (overlay cases and controls)
coords_neg <- cbind(matrix(rodents[rodents$Res=='NEG',]$Lon_Add_Fix), rodents[rodents$Res=='NEG',]$Lat_Add_Fix)
points(coords_neg, col=4, cex=0.2)
coords_pos <- cbind(matrix(rodents[rodents$Res=='POS',]$Lon_Add_Fix), rodents[rodents$Res=='POS',]$Lat_Add_Fix)
points(coords_pos, col=2, cex=0.2)

# export figure
par(mfrow=c(1,2))
plot(sigmaps$r_point_est)
plot(sigmaps$r_inds_95)


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

X_low <- load_x_ca(factor=5)
lodds_low_sp <- X_low %*% beta_ca_sp + w_ca_est - X_low %*% beta_co_sp - w_co_est
risk_low_sp <- calc_risk(lodds_low_sp)

# downscale the spatial poisson model via spline interpolation for speed
downscale <- function(w.est, caPr.disc, caPr){
  
  rw <- caPr.disc[[1]]
  rw[][!is.na(rw[])] <- w.est
  xy <- data.frame(xyFromCell(rw, 1:ncell(rw)))
  v <- getValues(rw)
  tps <- Tps(xy, v)
  p <- raster(caPr[[2]])
  p <- interpolate(p, tps)
  p <- mask(p, caPr[[1]])
  w.hat_ds <- p[][!is.na(p[])]
  return(list(p=p, w.hat_ds=w.hat_ds))
  
}

w_ca_ds <- downscale(w_ca_est, caPr.disc, caPr)$w.hat_ds
w_co_ds <- downscale(w_co_est, caPr.disc, caPr)$w.hat_ds

X_high <- load_x_ca()
lodds_high_sp <- X_high %*% beta_ca_sp + w_ca_ds - X_high %*% beta_co_sp - w_co_ds
risk_high_sp <- calc_risk(lodds_high_sp)

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
