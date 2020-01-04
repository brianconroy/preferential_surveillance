#############################
# Choose simulation parameter
# values
#############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=9)
n_values(caPr.disc[[1]])
plot(caPr.disc)


#### Simulate gaussian process
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 1, -0.5)
locs <- simBernoulliLocCov(caPr.disc, beta.samp, w=W, seed=56)
sum(locs$status)
hist(locs$w)
plot(caPr.disc[[1]])
points(locs$coords)

# none: alpha_+ = 0, alpha_- = 0
# how to justify "medium", "high" designations
# medium: 10, 25

# prevalences:
# low, medium, high
# 0.05, 0.15, 0.25

sampling <- 'none'
Alpha.case <- 0
Alpha.ctrl <- 0
prev <- 'high'
beta.case <- c(2.75, 0.75, 0.25)
beta.ctrl <- c(3.6, 1, 0.5)
cov.disc <- caPr.disc


calc_ps_contribution(cov.disc, locs, beta.case, Alpha.case, beta.ctrl, Alpha.ctrl, W)


#### Simulate counts given locations
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W, seed=42)
ctrl.data <- simConditionalGp2(cov.disc, locs, beta.ctrl, Alpha.ctrl, W, seed=40)
sum(case.data$y)/sum(case.data$y + ctrl.data$y)

sum(case.data$y)
sum(ctrl.data$y)
print(case.data$y)
print(ctrl.data$y)
print(round(case.data$y/(case.data$y + ctrl.data$y), 3))

params <- list(
  sampling=sampling,
  prevalence=prev,
  beta.case=beta.case,
  beta.ctrl=beta.ctrl,
  alpha.case=Alpha.case,
  alpha.ctrl=Alpha.ctrl,
  total.y.ca=sum(case.data$y),
  total.y.co=sum(ctrl.data$y),
  prev=round(sum(case.data$y)/sum(case.data$y + ctrl.data$y), 2),
  beta.samp=beta.samp,
  Theta=Theta,
  Phi=Phi,
  W=W
)

save_output(params, paste("v2_params_", prev, "_", sampling, ".json", sep=""))

data <- list(
  case.data=case.data,
  ctrl.data=ctrl.data,
  locs=locs
)

save_output(data, paste("v2_data_", prev, "_", sampling, ".json", sep=""))
