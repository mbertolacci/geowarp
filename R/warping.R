#' Construct a Bernstein Warping Design Matrix
#'
#' Creates a design matrix based on the Bernstein basis for use in GeoWarp's
#' axial warping units. The difference to the usual Bernstein basis is that the
#' first basis function is omitted so that input coordinate zero maps to warped
#' coordinate zero. The order determines the flexibility of the basis expansion.
#'
#' @param x Numeric vector representing the coordinate values; must be within
#' the given domain.
#' @param order An integer representing the order of the Bernstein polynomial.
#' @param domain A numeric vector of length 2 specifying the domain over which
#' the Bernstein polynomial is defined.
#'
#' @return A matrix where each column represents a Bernstein polynomial
#' evaluated at the data points `x`.
#'
#' @examples
#' # Create a Bernstein warping design matrix with 10 data points, of order 4,
#' # over domain [0, 5]
#' design_matrix <- bernstein_warping_design_matrix(
#'   x = seq(0, 5, length.out = 10),
#'   order = 4,
#'   domain = c(0, 5)
#' )
#' print(design_matrix)
#'
#' @export
bernstein_warping_design_matrix <- function(
  x,
  order,
  domain
) {
  x_prime <- (x - domain[1]) / (domain[2] - domain[1])
  stopifnot(all(x_prime >= 0 & x_prime <= 1))
  sapply(seq_len(order), function(k) {
    choose(order, k) * x_prime ^ k * (1 - x_prime) ^ (order - k)
  })
}

#' Warped and Unwarp Coordinates
#'
#' These function produces the warped coordinates for a GeoWarp model given the
#' original coordinates, or give the original coordinates given the warped
#' coordinates.
#'
#' @param fit A fitted GeoWarp object.
#' @param df A dataframe containing the original coordinate space data.
#' @param model A GeoWarp model object, defaulting to `fit$model`.
#' @param parameters Parameters for the GeoWarp model, defaulting to
#' `fit$parameters`.
#' @param x The unwarped coordinates as a matrix, by default equal to
#' `model_coordinates(df, fit$model)`.
#' @param warped_x Matrix containing the warped coordinates to be transformed
#' back.
#'
#' @return A matrix where each row represents a point in the warped coordinate
#' space.
#'
#' @seealso
#' \code{\link[geowarp]{geowarp_model}}
#' \code{\link[geowarp]{geowarp_optimize}}
#' \code{\link[geowarp]{model_coordinates}}
#'
#' @export
warped_coordinates <- function(
  fit,
  df,
  model = fit$model,
  parameters = fit$parameters,
  x = model_coordinates(df, model)
) {
  if (model$deviation_model$name == 'white') {
    colnames(x) <- NULL
    return(x)
  }
  is_vertical_only <- model$deviation_model$name == 'vertical_only'

  vertical_index <- length(model$horizontal_coordinates) + 1L
  output <- x
  if (!is_vertical_only) {
    for (i in seq_along(model$horizontal_coordinates)) {
      k <- model$deviation_model$axial_warping_unit_mapping[i]
      output[, i] <- (
        parameters$gamma_deviation_horizontal[k]
        * model$deviation_model$axial_warping_units[[k]]$scaling
        * output[, i]
      )
    }
  }

  vertical_warping <- if (is_vertical_only) {
    model$deviation_model$axial_warping_unit
  } else {
    tail(model$deviation_model$axial_warping_units, 1)[[1]]
  }
  if (vertical_warping$name == 'linear_awu') {
    output[, vertical_index] <- (
      parameters$gamma_deviation_vertical[1]
      * vertical_warping$scaling
      * output[, vertical_index]
    )
  } else {
    X_deviation_warping <- vertical_warping$scaling * bernstein_warping_design_matrix(
      output[, vertical_index],
      vertical_warping$order,
      model$vertical_domain
    )
    output[, vertical_index] <- (
      X_deviation_warping
      %*% cumsum(parameters$gamma_deviation_vertical)
    )
  }

  if (!is.null(model$deviation_model$geometric_warping_unit)) {
    output <- output %*% parameters$L_deviation
  }

  colnames(output) <- NULL

  output
}

#' @describeIn warped_coordinates Calculated unwarped coordinates.
#' @export
unwarped_coordinates <- function(
  fit,
  warped_x,
  model = fit$model,
  parameters = fit$parameters
) {
  if (model$deviation_model$name == 'white') {
    return(warped_x)
  }

  is_vertical_only <- model$deviation_model$name == 'vertical_only'
  vertical_index <- length(model$horizontal_coordinates) + 1L
  output <- warped_x

  if (!is.null(model$deviation_model$geometric_warping_unit)) {
    output <- output %*% solve(parameters$L_deviation)
  }

  if (!is_vertical_only) {
    for (i in seq_along(model$horizontal_coordinates)) {
      k <- model$deviation_model$axial_warping_unit_mapping[i]
      output[, i] <- (
        output[, i]
        / parameters$gamma_deviation_horizontal[k]
        / model$deviation_model$axial_warping_units[[k]]$scaling
      )
    }
  }

  vertical_warping <- if (is_vertical_only) {
    model$deviation_model$axial_warping_unit
  } else {
    tail(model$deviation_model$axial_warping_units, 1)[[1]]
  }
  if (vertical_warping$name == 'linear_awu') {
    output[, vertical_index] <- (
      output[, vertical_index]
      / parameters$gamma_deviation_vertical[1]
      / vertical_warping$scaling
    )
  } else {
    intervals <- matrix(0, nrow = nrow(output), ncol = 2)
    intervals[, 1] <- model$vertical_domain[1]
    intervals[, 2] <- model$vertical_domain[2]
    for (iteration in seq_len(1000)) {
      x_mid <- rowMeans(intervals)
      y_mid <- vertical_warping$scaling * bernstein_warping_design_matrix(
        x_mid,
        vertical_warping$order,
        model$vertical_domain
      ) %*% cumsum(parameters$gamma_deviation_vertical)

      is_negative <- y_mid - output[, vertical_index] < 0
      intervals[is_negative, 1] <- x_mid[is_negative]
      intervals[!is_negative, 2] <- x_mid[!is_negative]

      if (all(intervals[, 2] - intervals[, 1] < .Machine$double.eps ^ 0.5)) {
        break
      }
    }
    output[, vertical_index] <- rowMeans(intervals)
  }

  output
}
