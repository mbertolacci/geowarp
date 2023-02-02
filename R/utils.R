.chol_solve <- function(R, b) {
  backsolve(R, backsolve(R, b, transpose = TRUE))
}

#' @export
gamma_quantile_prior <- function(
  q_lower, q_upper,
  p_lower = 0.05, p_upper = 0.95,
  interval = c(0.0001, 1000)
) {
  stopifnot(q_lower < q_upper)
  stopifnot(p_lower < p_upper)

  # Find the shape parameter that gives the appropriate ratio between the
  # quantiles
  ratio <- q_lower / q_upper
  shape <- uniroot(function(shape) {
    theoretical <- qgamma(c(p_lower, p_upper), shape = shape, rate = 1)
    theoretical[1] / theoretical[2] - ratio
  }, interval, tol = sqrt(.Machine$double.eps))$root

  # Find the rate parameter that gives the correct quantiles
  rate <- qgamma(p_upper, shape = shape, rate = 1) / q_upper

  list(shape = shape, rate = rate)
}

#' @export
inverse_gamma_quantile_prior <- function(
  q_lower, q_upper, ...
) {
  gamma_quantile_prior(1 / q_upper, 1 / q_lower, ...)
}
