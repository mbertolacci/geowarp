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
  warping = pcpt_linear_warping(),
  variance_model = pcpt_variance_model()
) {
  list(
    name = 'vertical_only',
    warping = warping,
    variance_model = variance_model
  )
}

#' @export
pcpt_isotropic_deviation_model <- function(
  warpings = list(
    pcpt_linear_warping(),
    pcpt_linear_warping(),
    pcpt_linear_warping()
  ),
  variance_model = pcpt_variance_model()
) {
  .check_horizontal_warping(warpings)
  list(
    name = 'isotropic',
    warpings = warpings,
    variance_model = variance_model
  )
}

#' @export
pcpt_anisotropic_deviation_model <- function(
  warpings = list(
    pcpt_linear_warping(),
    pcpt_linear_warping(),
    pcpt_linear_warping()
  ),
  anisotropy_shape = 1,
  variance_model = pcpt_variance_model()
) {
  .check_horizontal_warping(warpings)
  list(
    name = 'anisotropic',
    warpings = warpings,
    anisotropy_shape = anisotropy_shape,
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
pcpt_linear_warping <- function(
  scaling = 1,
  prior = list(shape = 1, rate = 0)
) {
  list(name = 'linear', scaling = scaling, prior = prior)
}

#' @export
pcpt_bernstein_warping <- function(
  order = 20,
  domain,
  prior = list(shape = 1, rate = 0)
) {
  stopifnot(domain[2] > domain[1])
  list(
    name = 'bernstein_warping',
    order = order,
    scaling = 1,
    domain = domain,
    prior = prior
  )
}


.check_horizontal_warping <- function(warpings) {
  warping_names <- sapply(warpings, getElement, 'name')
  if (!all(head(warping_names, -1) == 'linear')) {
    stop('non-linear warping currently not supported for horizontal dimensions')
  }
}


#' @export
pcpt_length_scale_prior <- function(
  q_lower, q_upper,
  p_lower = 0.05, p_upper = 0.95
) {
  gamma_quantile_prior(
    1 / q_upper, 1 / q_lower,
    p_lower, p_upper
  )
}
