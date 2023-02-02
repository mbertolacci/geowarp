# #' @export
# arial_warping_unit_design_matrix <- function(
#   x,
#   n_basis_functions,
#   gradient,
#   domain
# ) {
#   locations <- seq(domain[1], domain[2], length.out = n_basis_functions)
#   do.call(cbind, lapply(locations, function(location_i) {
#     1 / (1 + exp(-gradient * (x - location_i)))
#   }))
# }

#' @export
bernstein_warping_design_matrix <- function(
  x,
  order,
  domain
) {
  x_prime <- (x - domain[1]) / (domain[2] - domain[1])
  stopifnot(min(x_prime) >= 0 && max(x_prime) <= 1)
  # sapply(seq_len(order), function(k) {
  #   pracma::bernsteinb(k, order, x_prime)
  # })
  sapply(0 : order, function(k) {
    pracma::bernsteinb(k, order, x_prime)
  })
}

#' @export
warped_coordinates <- function(
  fit,
  df,
  x = as.matrix(
    df[, c(fit$model$horizontal_coordinates, fit$model$vertical_coordinate)]
  )
) {
  model <- fit$model
  if (model$deviation_model$name == 'white') {
    return(x)
  }
  is_vertical_only <- model$deviation_model$name == 'vertical_only'

  vertical_index <- length(model$horizontal_coordinates) + 1L
  output <- x
  if (!is_vertical_only) {
    for (i in seq_along(model$horizontal_coordinates)) {
      output[, i] <- (
        fit$parameters$gamma_deviation_horizontal[i]
        * model$deviation_model$warpings[[i]]$scaling
        * output[, i]
      )
    }
  }

  vertical_warping <- if (is_vertical_only) {
    model$deviation_model$warping
  } else {
    tail(model$deviation_model$warpings, 1)[[1]]
  }
  if (vertical_warping$name == 'linear') {
    output[, vertical_index] <- (
      fit$parameters$gamma_deviation_vertical[1]
      * vertical_warping$scaling
      * output[, vertical_index]
    )
  } else {
    X_deviation_warping <- bernstein_warping_design_matrix(
      vertical_warping$scaling * output[, vertical_index],
      vertical_warping$order,
      vertical_warping$domain
    )
    output[, vertical_index] <- X_deviation_warping %*% cumsum(fit$parameters$gamma_deviation_vertical)
  }

  if (model$deviation_model$name == 'anisotropic') {
    output <- output %*% fit$parameters$L_deviation
  }

  output
}
