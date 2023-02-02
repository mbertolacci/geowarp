.cor_isotropic_exp <- function(x) {
  exp(-1.0 * fields::rdist(x))
}

.cor_anisotropic_exp <- function(x, ell, L) {
  exp(-1.0 * fields::rdist(
    x %*% (diag(1 / ell) %*% L)
  ))
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

# #' @export
# pcpt_covariance_matrix <- function(
#   fit,
#   df = fit$input_df,
#   model = fit$model,
#   parameters
# ) {
#   .covariance_matrix_internal(
#     .to_stan_data(df, model),
#     model,
#     parameters
#   )
# }

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
  nugget
) {
  sigma_deviation <- exp(
    0.5 * as.vector(
      X_deviation_fixed %*% fit$parameters$eta_deviation
      + X_deviation_random %*% fit$parameters$zeta_deviation
    )
  )
  if (missing(nugget)) {
    nugget <- rep(TRUE, nrow(x))
  }

  x_warped <- warped_coordinates(fit, x = x)
  deviation_cor <- if (fit$model$deviation_model$name == 'vertical_only') {
    .cor_vertical_only_exp(x_warped)
  } else {
    .cor_isotropic_exp(x_warped)
  }

  output <- .quad_form_diag(deviation_cor, sigma_deviation)
  diag(output) <- diag(output) + nugget * fit$parameters$sigma_squared_nugget

  output
}
