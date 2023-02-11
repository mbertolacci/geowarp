.cor_isotropic_exp <- function(x) {
  exp(-1.0 * fields::rdist(x))
}

.cor_vertical_only_exp <- function(x) {
  output <- matrix(0, nrow = nrow(x), ncol = nrow(x))

  stopifnot(ncol(x) <= 3)
  if (ncol(x) == 3) {
    location_coords <- paste0(x[, 1], x[, 2])
  } else {
    location_coords <- x[, 1]
  }
  for (location_coord_i in unique(location_coords)) {
    indices <- location_coords == location_coord_i
    output[indices, indices] <- .cor_isotropic_exp(x[indices, ncol(x)])
  }
  output
}

.quad_form_diag <- function(Y, x) {
  output <- Y
  output[] <- x * Y * rep(x, each = ncol(Y))
  output
}

.covariance_matrix_internal <- function(
  x,
  X_deviation_fixed,
  X_deviation_random,
  fit,
  nugget,
  parameters = fit$parameters,
  model = fit$model
) {
  sigma_deviation <- exp(
    0.5 * as.vector(
      X_deviation_fixed %*% parameters$eta_deviation
      + X_deviation_random %*% parameters$zeta_deviation
    )
  )
  if (missing(nugget)) {
    nugget <- rep(TRUE, nrow(x))
  }

  x_warped <- warped_coordinates(x = x, model = model, parameters = parameters)
  deviation_cor <- if (model$deviation_model$name == 'vertical_only') {
    .cor_vertical_only_exp(x_warped)
  } else {
    .cor_isotropic_exp(x_warped)
  }

  output <- .quad_form_diag(deviation_cor, sigma_deviation)
  diag(output) <- diag(output) + nugget * parameters$sigma_squared_nugget

  output
}
