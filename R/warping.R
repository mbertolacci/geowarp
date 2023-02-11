#' @export
bernstein_warping_design_matrix <- function(
  x,
  order,
  domain
) {
  x_prime <- (x - domain[1]) / (domain[2] - domain[1])
  stopifnot(min(x_prime) >= 0 && max(x_prime) <= 1)
  sapply(seq_len(order), function(k) {
    pracma::bernsteinb(k, order, x_prime)
  })
}

#' @export
warped_coordinates <- function(
  fit,
  df,
  model = fit$model,
  parameters = fit$parameters,
  x = as.matrix(
    df[, c(model$horizontal_coordinates, model$vertical_coordinate)]
  )
) {
  if (model$deviation_model$name == 'white') {
    return(x)
  }
  is_vertical_only <- model$deviation_model$name == 'vertical_only'

  vertical_index <- length(model$horizontal_coordinates) + 1L
  output <- x
  if (!is_vertical_only) {
    for (i in seq_along(model$horizontal_coordinates)) {
      output[, i] <- (
        parameters$gamma_deviation_horizontal[i]
        * model$deviation_model$axial_warping_units[[i]]$scaling
        * output[, i]
      )
    }
  }

  vertical_warping <- if (is_vertical_only) {
    model$deviation_model$axial_warping_unit
  } else {
    tail(model$deviation_model$axial_warping_units, 1)[[1]]
  }
  if (vertical_warping$name == 'linear') {
    output[, vertical_index] <- (
      parameters$gamma_deviation_vertical[1]
      * vertical_warping$scaling
      * output[, vertical_index]
    )
  } else {
    X_deviation_warping <- bernstein_warping_design_matrix(
      vertical_warping$scaling * output[, vertical_index],
      vertical_warping$order,
      vertical_warping$domain
    )
    output[, vertical_index] <- (
      X_deviation_warping
      %*% cumsum(parameters$gamma_deviation_vertical)
    )
  }

  if (!is.null(model$deviation_model$rotation_unit)) {
    output <- output %*% parameters$L_deviation
  }

  output
}
