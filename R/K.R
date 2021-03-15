

#' K
#' 
#' Returns the kinetic energy of a momentum vector. Used as an auxiliary
#' function in Hamiltonian Monte Carlo samplers.
#'
#' @param p numeric. vector of momenta.
#'
#' @return numeric. kinetic energy.
#' @export
#'
#' @examples
#' 
#' momentum <- rnorm(20)
#' K(momentum)
#' 
K <- function(p){
  # kinetic energy
  return(t(p) %*% p/2)
}
