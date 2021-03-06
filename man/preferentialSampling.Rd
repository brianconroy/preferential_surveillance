% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preferentialSampling.R
\name{preferentialSampling}
\alias{preferentialSampling}
\title{preferentialSampling}
\usage{
preferentialSampling(
  data,
  d,
  n.sample,
  burnin,
  L_w,
  L_ca,
  L_co,
  L_a_ca,
  L_a_co,
  proposal.sd.theta = 0.3,
  m_aca = 2000,
  m_aco = 2000,
  m_ca = 700,
  m_co = 700,
  m_w = 700,
  target_aca = 0.75,
  target_aco = 0.75,
  target_ca = 0.75,
  target_co = 0.75,
  target_w = 0.75,
  target_loc = 0.75,
  self_tune_w = TRUE,
  self_tune_aca = TRUE,
  self_tune_aco = TRUE,
  self_tune_ca = TRUE,
  self_tune_co = TRUE,
  self_tune_loc = TRUE,
  delta_w = NULL,
  delta_aca = NULL,
  delta_aco = NULL,
  delta_ca = NULL,
  delta_co = NULL,
  delta_loc = NULL,
  beta_ca_initial = NULL,
  beta_co_initial = NULL,
  alpha_ca_initial = NULL,
  alpha_co_initial = NULL,
  beta_loc_initial = NULL,
  theta_initial = NULL,
  phi_initial = NULL,
  w_initial = NULL,
  prior_phi,
  prior_theta,
  prior_alpha_ca_var,
  prior_alpha_co_var
)
}
\arguments{
\item{data}{List. Data input containing case, control counts, covariates.}

\item{d}{Matrix. Distance matrix for grid cells in study region.}

\item{n.sample}{Numeric. Number of MCMC samples to generate.}

\item{burnin}{Numeric. Number of MCMC samples to discard as burnin.}

\item{L_w}{Numeric. HMC simulation length parameter for spatial random effects.}

\item{L_ca}{Numeric. HMC simulation length parameter for case covariates.}

\item{L_co}{Numeric. HMC simulation length parameter for control covariates.}

\item{L_a_ca}{Numeric. HMC simulation length parameter for case preferential sampling parameter.}

\item{L_a_co}{Numeric. HMC simulation length parameter for control preferential sampling parameter.}

\item{proposal.sd.theta}{Numeric. Standard deviation of proposal distribution for spatial range parameter.}

\item{m_aca}{Numeric. Number of samples to apply self tuning for case preferential sampling parameter.}

\item{m_aco}{Numeric. Number of samples to apply self tuning for control preferential sampling parameter.}

\item{m_ca}{Numeric. Number of samples to apply self tuning for case covariates.}

\item{m_co}{Numeric. Number of samples to apply self tuning for control covariates.}

\item{m_w}{Numeric. Number of samples to apply self tuning for spatial random effects.}

\item{target_aca}{Numeric. Target acceptance rate for case preferential sampling parameter.}

\item{target_aco}{Numeric. Target acceptance rate for control preferential sampling parameter.}

\item{target_ca}{Numeric. Target acceptance rate for case covariates.}

\item{target_co}{Numeric. Target acceptance rate for control covariates.}

\item{target_w}{Numeric. Target acceptance rate for spatial random effects.}

\item{target_loc}{Numeric. Target acceptance rate for locational covariates.}

\item{self_tune_w}{Logical. Whether to apply self tuning for spatial random effects.}

\item{self_tune_aca}{Logical. Whether to apply self tuning for case preferential sampling paramter.}

\item{self_tune_aco}{Logical. Whether to apply self tuning for control preferential sampling paramter.}

\item{self_tune_ca}{Logical. Whether to apply self tuning for case covariates.}

\item{self_tune_co}{Logical. Whether to apply self tuning for control covariates.}

\item{self_tune_loc}{Logical. Whether to apply self tuning for locational covariates.}

\item{delta_w}{Numeric. Required if self_tune_w is FALSE. HMC step size for spatial random effects.}

\item{delta_aca}{Numeric. Required if self_tune_w is FALSE. HMC step size for case preferential sampling parameter.}

\item{delta_aco}{Numeric. Required if self_tune_w is FALSE. HMC step size for control preferential sampling parameter.}

\item{delta_ca}{Numeric. Required if self_tune_w is FALSE. HMC step size for case covariates.}

\item{delta_co}{Numeric. Required if self_tune_w is FALSE. HMC step size for control covariates.}

\item{delta_loc}{Numeric. Required if self_tune_w is FALSE. HMC step size for locational covariates.}

\item{beta_ca_initial}{Numeric. Initial MCMC value for case covariate parameter.}

\item{beta_co_initial}{Numeric. Initial MCMC value for control covariate parameter.}

\item{alpha_ca_initial}{Numeric. Initial MCMC value for case preferential sampling parameter.}

\item{alpha_co_initial}{Numeric. Initial MCMC value for control preferential sampling parameter.}

\item{beta_loc_initial}{Numeric. Initial MCMC value for locational covariate parameter.}

\item{theta_initial}{Numeric. Initial MCMC value for spatial range parameter.}

\item{phi_initial}{Numeric. Initial MCMC value for spatial marginal variance parameter.}

\item{w_initial}{Numeric. Initial MCMC value for spatial random effects.}

\item{prior_phi}{List. Shape and scale values for prior distribution (Inverse Gamma) of marginal variance.}

\item{prior_theta}{List. Shape and scale values for prior distribution (Gamma) of spatial range.}

\item{prior_alpha_ca_var}{List. Prior (Independent Normal) variance of case covariates.}

\item{prior_alpha_co_var}{List. Prior (Independent Normal) variance of control covariates.}
}
\value{
List containing posterior samples and associated tuning values.
}
\description{
preferentialSampling
}
