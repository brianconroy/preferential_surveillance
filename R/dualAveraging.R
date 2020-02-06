

initialize_tuning <- function(m, target){
  return(list(
    M_adapt=m,
    delta_curr=0.05,
    delta_bar=1,
    delta_tar=target,
    mu=log(10*0.05),
    HbarM=0,
    gamma=0.05,
    t0=10,
    kappa=0.75,
    reject_streak=0
  ))
}


update_tuning <- function(tuning, a, i, accepted){
  
  if (accepted){
    tuning$reject_streak <- 0
  } else{
    tuning$reject_streak <- tuning$reject_streak + 1
  }
  
  if (i <= tuning$M_adapt){
    
    delta_tar <- tuning$delta_tar
    delta_bar <- tuning$delta_bar
    HbarM <- tuning$HbarM
    gamma <- tuning$gamma
    mu <- tuning$mu
    kappa <- tuning$kappa
    t0 <- tuning$t0
    
    HbarM <- (1 - 1/(i + t0)) * HbarM + (1/(i + t0)) * (delta_tar - a)
    log_delta_curr <- mu - (sqrt(i)/gamma) * HbarM
    log_delta_bar <- i^{-kappa} * log_delta_curr + (1 - i^{-kappa}) * log(delta_bar)
    delta_curr <- exp(log_delta_curr)
    delta_bar <- exp(log_delta_bar)
    
    tuning$delta_bar <- delta_bar
    tuning$HbarM <- HbarM
    tuning$delta_curr <- delta_curr
    
  } else{
    
    tuning$delta_curr <- tuning$delta_bar
    if (tuning$reject_streak > 1000) {
      tuning$delta_curr <- 0.90 * tuning$delta_curr
      tuning$delta_bar <- 0.90 * tuning$delta_bar
    }
    
  }
  
  return(tuning)
  
}
