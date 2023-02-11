#' @export
pcpt_model <- function(
  variable,
  horizontal_coordinates,
  vertical_coordinate,
  mean_model = pcpt_mean_model(),
  deviation_model,
  nugget_prior = list(shape = 1, rate = 0)
) {
  output <- structure(list(
    variable = variable,
    horizontal_coordinates = horizontal_coordinates,
    vertical_coordinate = vertical_coordinate,
    mean_model = mean_model,
    deviation_model = deviation_model,
    nugget_prior = nugget_prior
  ), class = 'pcpt_model')
  output
}

#' @export
print.pcpt_model <- function(x, ...) {
  cat(sprintf(
    'PCPT model with %d horizontal dimensions and a deviation model that is %s\n',
    length(x$horizontal_coordinates),
    x$deviation_model$name
  ))
  invisible(x)
}

#' @export
pcpt_mean_model <- function(
  fixed_formula = ~ 1,
  fixed_effect_mean = 0,
  fixed_effect_precision = 1 / 100,
  vertical_basis_functions = FALSE,
  vertical_basis_function_domain = c(0, 40),
  vertical_basis_function_delta = 1,
  vertical_basis_function_effect_precision = 1 / 100,
  vertical_basis_function_variance_prior = list(scale = 0)
) {
  list(
    fixed_formula = fixed_formula,
    fixed_effect_mean = fixed_effect_mean,
    fixed_effect_precision = fixed_effect_precision,
    vertical_basis_functions = vertical_basis_functions,
    vertical_basis_function_domain = vertical_basis_function_domain,
    vertical_basis_function_delta = vertical_basis_function_delta,
    vertical_basis_function_effect_precision = vertical_basis_function_effect_precision,
    vertical_basis_function_variance_prior = vertical_basis_function_variance_prior
  )
}

#' @export
pcpt_white_deviation_model <- function(
  variance_model = pcpt_variance_model()
) {
  list(
    name = 'white',
    variance_model = variance_model
  )
}

#' @export
pcpt_vertical_only_deviation_model <- function(
  axial_warping_unit = pcpt_linear_awu(),
  variance_model = pcpt_variance_model()
) {
  list(
    name = 'vertical_only',
    axial_warping_unit = axial_warping_unit,
    variance_model = variance_model
  )
}

#' @export
pcpt_full_deviation_model <- function(
  axial_warping_units = list(
    pcpt_linear_awu(),
    pcpt_linear_awu(),
    pcpt_linear_awu()
  ),
  rotation_unit = NULL,
  variance_model = pcpt_variance_model()
) {
  .check_horizontal_warping(axial_warping_units)

  list(
    name = 'full',
    axial_warping_units = axial_warping_units,
    rotation_unit = rotation_unit,
    variance_model = variance_model
  )
}

#' @export
pcpt_variance_model <- function(
  fixed_formula = ~ 1,
  fixed_effect_mean = 0,
  fixed_effect_precision = 1 / 100,
  vertical_basis_functions = FALSE,
  vertical_basis_function_domain = c(0, 40),
  vertical_basis_function_delta = 1,
  vertical_basis_function_length_scale_prior = list(scale = 0),
  vertical_basis_function_variance_prior = list(scale = 0)
) {
  if (vertical_basis_functions) {
    stopifnot(vertical_basis_function_domain[2] > vertical_basis_function_domain[1])
  }
  list(
    fixed_formula = fixed_formula,
    fixed_effect_mean = fixed_effect_mean,
    fixed_effect_precision = fixed_effect_precision,
    vertical_basis_functions = vertical_basis_functions,
    vertical_basis_function_domain = vertical_basis_function_domain,
    vertical_basis_function_delta = vertical_basis_function_delta,
    vertical_basis_function_length_scale_prior = vertical_basis_function_length_scale_prior,
    vertical_basis_function_variance_prior = vertical_basis_function_variance_prior
  )
}

#' @export
pcpt_linear_awu <- function(
  scaling = 1,
  prior = list(shape = 1, rate = 0)
) {
  list(name = 'linear_awu', scaling = scaling, prior = prior)
}

#' @export
pcpt_bernstein_awu <- function(
  order = 20,
  domain,
  prior = list(shape = 1, rate = 0)
) {
  stopifnot(domain[2] > domain[1])
  list(
    name = 'bernstein_awu',
    order = order,
    scaling = 1,
    domain = domain,
    prior = prior
  )
}

#' @export
pcpt_rotation_warping_unit <- function(prior_shape = 1) {
  list(
    name = 'rotation_warping_unit',
    prior_shape = prior_shape
  )
}

.check_horizontal_warping <- function(warpings) {
  warping_names <- sapply(warpings, getElement, 'name')
  if (!all(head(warping_names, -1) == 'linear_awu')) {
    stop('non-linear warping currently not supported for horizontal dimensions')
  }
}
