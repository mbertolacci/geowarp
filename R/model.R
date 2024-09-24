#' Construct a GeoWarp Model Object
#'
#' Creates a GeoWarp model object that specifies the structure and priors of
#' the model. The created object is intended to be used as input for
#' \code{\link{geowarp_optimise}}.
#'
#' @param variable Character vector naming the target variable.
#' @param horizontal_coordinates Vector containing names of horizontal
#' coordinates.
#' @param vertical_coordinate Vector containing name of vertical coordinate.
#' @param horizontal_domains Domain of horizontal coordinates, or a list if
#' there is more than one. A domain is a numeric vector of length two.
#' @param vertical_domain Domain of the vertical coordinate, a vector of length
#' two.
#' @param mean_model GeoWarp mean model object. Defaults to the model given in
#' the GeoWarp paper, with an intercept and a trend based on the vertical
#' coordinate. See \code{\link{geowarp_mean_model}} for details.
#' @param deviation_model GeoWarp deviation model object. Defaults match the
#' configuration given in the GeoWarp paper.
#' @param nugget_prior Prior distribution for the nugget effect given as a list
#' with `shape` and `rate` parameters.
#'
#' @return A `geowarp_model` object containing the specifications for the
#' GeoWarp model, ready to be passed to \code{\link{geowarp_optimise}}.
#'
#' @examples
#' # Creating a simple GeoWarp model object for a 2-D problem
#' model <- geowarp_model(
#'   variable = 'log_q_c',
#'   horizontal_coordinates = 'horizontal',
#'   vertical_coordinate = 'depth',
#'   horizontal_domains = c(0, 150),
#'   vertical_domain = c(0, 43)
#' )
#' # Creating a simple GeoWarp model object for a 3-D problem
#' model <- geowarp_model(
#'   variable = 'log_q_c',
#'   horizontal_coordinates = c('easting', 'northing'),
#'   vertical_coordinate = 'depth',
#'   horizontal_domains = list(c(0, 150), c(0, 200)),
#'   vertical_domain = c(0, 43)
#' )
#' # Overriding some of the parameters, namely the spacing of the splines used
#' # in both the mean model and the depth-varying variances
#' model <- geowarp_model(
#'   variable = 'log_q_c',
#'   horizontal_coordinates = c('easting', 'northing'),
#'   vertical_coordinate = 'depth',
#'   horizontal_domains = list(c(0, 150), c(0, 200)),
#'   vertical_domain = c(0, 43),
#'   mean_model = geowarp_mean_model(
#'     fixed_formula = ~ depth,
#'     # Set the spacing for knots in the mean function to 1
#'     vertical_basis_function_delta = 1
#'   ),
#'   deviation_model = geowarp_deviation_model(
#'     axial_warping_units = list(
#'       geowarp_linear_awu(),
#'       geowarp_linear_awu(),
#'       geowarp_bernstein_awu()
#'     ),
#'     variance_model = geowarp_variance_model(
#'       # Set knot spacing to 5
#'       vertical_basis_function_delta = 5
#'     )
#'   )
#' )
#'
#' @seealso
#' \code{\link{geowarp_mean_model}}
#' \code{\link{geowarp_deviation_model}}
#' \code{\link{geowarp_optimise}}
#' @export
geowarp_model <- function(
  variable,
  horizontal_coordinates,
  vertical_coordinate,
  horizontal_domains,
  vertical_domain,
  mean_model = geowarp_mean_model(
    fixed_formula = vertical_coordinate
  ),
  deviation_model = geowarp_deviation_model(
    axial_warping_units = c(
      rep(list(geowarp_linear_awu()), length(horizontal_coordinates)),
      list(geowarp_bernstein_awu())
    )
  ),
  nugget_prior = inverse_gamma_quantile_prior(0.1, 1)
) {
  if (!is.list(horizontal_domains)) {
    horizontal_domains <- list(horizontal_domains)
  }

  stopifnot(length(horizontal_coordinates) == length(horizontal_domains))
  for (horizontal_domain in horizontal_domains) {
    stopifnot(horizontal_domain[2] > horizontal_domain[1])
  }
  stopifnot(vertical_domain[2] > vertical_domain[1])

  structure(list(
    variable = variable,
    horizontal_coordinates = horizontal_coordinates,
    horizontal_domains = horizontal_domains,
    vertical_coordinate = vertical_coordinate,
    vertical_domain = vertical_domain,
    mean_model = mean_model,
    deviation_model = deviation_model,
    nugget_prior = nugget_prior
  ), class = 'geowarp_model')
}

#' Create a GeoWarp Mean Model Object
#'
#' Constructs a GeoWarp mean model object that specifies the fixed effects
#' formula, the mean and precision of the fixed effects, and parameters for
#' vertical basis functions. This object is intended to be used as an argument
#' for \code{\link{geowarp_model}}.
#'
#' @param fixed_formula A formula or character vector that specifies the fixed
#' effects. If a character vector is provided, it will be assumed to refer to
#' a single variable and converted to a formula.
#' @param fixed_effect_mean A numeric value representing the mean of the fixed
#' effect. Defaults to 0.
#' @param fixed_effect_precision A numeric value representing the precision of
#' the fixed effect. Defaults to 1 / 100.
#' @param vertical_basis_functions Logical value indicating whether to include
#' vertical basis functions. Defaults to TRUE.
#' @param vertical_basis_function_delta A numeric value that specifies the
#' spacing of the vertical basis functions. Defaults to 0.1.
#' @param vertical_basis_function_effect_precision A numeric value representing
#' the precision of the effect of the vertical basis functions. Defaults to
#' 1 / 100.
#' @param vertical_basis_function_variance_prior Prior distribution for the
#' variance of the vertical basis functions, given as a list with `shape` and
#' `rate` parameters. Defaults to an inverse gamma prior with 5% and 95%
#' percentiles of `1e-3^2` and `10^2` respectively
#' @param vertical_basis_function_boundary_knots The number of knots to place
#' before after after the domain in order to address boundary effects. Defaults
#' to three.
#'
#' @return A GeoWarp mean model, ready to be passed as an argument to
#' \code{\link{geowarp_model}}.
#'
#' @examples
#' # Creating a GeoWarp mean model object with custom settings
#' mean_model <- geowarp_mean_model(
#'   fixed_formula = ~ depth,
#'   fixed_effect_mean = 1,
#'   fixed_effect_precision = 1,
#'   vertical_basis_function_delta = 0.2
#' )
#'
#' @seealso
#' \code{\link{geowarp_model}}, \code{\link{geowarp_deviation_model}}
#' @export
geowarp_mean_model <- function(
  fixed_formula,
  fixed_effect_mean = 0,
  fixed_effect_precision = 1 / 100,
  vertical_basis_functions = TRUE,
  vertical_basis_function_delta = 0.1,
  vertical_basis_function_effect_precision = 1 / 100,
  vertical_basis_function_variance_prior = inverse_gamma_quantile_prior(
    1e-3 ^ 2,
    10 ^ 2
  ),
  vertical_basis_function_boundary_knots = 3L
) {
  if (is.character(fixed_formula)) {
    fixed_formula <- as.formula(paste('~', fixed_formula))
  }
  list(
    fixed_formula = fixed_formula,
    fixed_effect_mean = fixed_effect_mean,
    fixed_effect_precision = fixed_effect_precision,
    vertical_basis_functions = vertical_basis_functions,
    vertical_basis_function_delta = vertical_basis_function_delta,
    vertical_basis_function_effect_precision = vertical_basis_function_effect_precision,
    vertical_basis_function_variance_prior = vertical_basis_function_variance_prior,
    vertical_basis_function_boundary_knots = vertical_basis_function_boundary_knots
  )
}

#' Create a GeoWarp Deviation Model Object
#'
#' Constructs a GeoWarp deviation model object specifying the correlation
#' structure, the axial warping units, and the variance model. This object is
#' intended to be used as an argument for \code{\link{geowarp_model}}.
#'
#' @param covariance_function A character vector specifying the type of
#' covariance function to use. Options are 'matern15' and 'exponential'.
#' Defaults to 'matern15'.
#' @param axial_warping_units A list of axial warping unit objects that specify
#' the warping functions for all the coordinates (for
#' `geowarp_deviation_model`). One of \code{\link{geowarp_linear_awu}}, or
#' \code{\link{geowarp_bernstein_awu}}. At the moment, horizontal coordinates
#' must only have linear warpings.
#' @param axial_warping_unit_mapping Specifies the mapping between the axial
#' warping units and the coordinates. Can be used to share warping units between
#' coordinates.
#' @param geometric_warping_unit A geometric warping unit object. Defaults to
#' `geowarp_geometric_warping_unit()`.
#' @param variance_model A variance model object. Defaults to
#' `geowarp_variance_model()`.
#'
#' @return A list containing the specifications for the GeoWarp deviation model,
#' ready to be passed as an argument to \code{\link{geowarp_model}}.
#'
#' @examples
#' # Creating a GeoWarp deviation model object with custom settings
#' deviation_model <- geowarp_deviation_model(
#'   covariance_function = 'exponential',
#'   axial_warping_units = list(geowarp_linear_awu(), geowarp_bernstein_awu())
#' )
#'
#' @seealso
#' \code{\link{geowarp_model}}, \code{\link{geowarp_mean_model}},
#' \code{\link{geowarp_linear_awu}}, \code{\link{geowarp_bernstein_awu}},
#' \code{\link{geowarp_geometric_warping_unit}},
#' \code{\link{geowarp_variance_model}}
#' @export
geowarp_deviation_model <- function(
  covariance_function = c('matern15', 'exponential'),
  axial_warping_units,
  axial_warping_unit_mapping = seq_len(length(axial_warping_units)),
  geometric_warping_unit = geowarp_geometric_warping_unit(),
  variance_model = geowarp_variance_model()
) {
  covariance_function <- match.arg(covariance_function)

  .check_horizontal_warping(axial_warping_units)

  list(
    name = 'full',
    covariance_function = covariance_function,
    axial_warping_units = axial_warping_units,
    axial_warping_unit_mapping = axial_warping_unit_mapping,
    geometric_warping_unit = geometric_warping_unit,
    variance_model = variance_model
  )
}

#' @describeIn geowarp_deviation_model A deviation model that assumes no
#' correlation (that is, i.i.d. Gaussian)
#' @export
geowarp_white_deviation_model <- function(
  variance_model = geowarp_variance_model()
) {
  list(
    name = 'white',
    variance_model = variance_model
  )
}

#' @describeIn geowarp_deviation_model A deviation model that assumes no
#' correlation in the horizontal directions.
#' @export
geowarp_vertical_only_deviation_model <- function(
  covariance_function = c('matern15', 'exponential'),
  axial_warping_unit = geowarp_bernstein_awu(),
  variance_model = geowarp_variance_model()
) {
  covariance_function <- match.arg(covariance_function)
  list(
    name = 'vertical_only',
    covariance_function = covariance_function,
    axial_warping_unit = axial_warping_unit,
    variance_model = variance_model
  )
}

#' Create a GeoWarp Variance Model Object
#'
#' Constructs a GeoWarp variance model object that specifies the model for
#' the log variance of the deviation process. The log variance has a fixed
#' effect, and optional vertical basis functions. This object is intended to be
#' used as an argument for \code{\link{geowarp_deviation_model}} and
#' \code{\link{geowarp_model}}.
#'
#' @param fixed_formula A formula specifying the fixed effect. Defaults to an
#' intercept-only model (`~ 1`).
#' @param fixed_effect_mean Numeric value for the prior mean of the fixed
#' effect.
#' @param fixed_effect_precision Numeric value for the prior precision of the
#' fixed effect.
#' @param vertical_basis_functions Logical indicating whether to include
#' vertical basis functions. Defaults to `TRUE`.
#' @param vertical_basis_function_delta Numeric value specifying the spacing
#' between knots for the vertical basis functions. Defaults to 1.
#' @param vertical_basis_function_length_scale_prior List containing the prior
#' for the length scale parameter for the basis function coefficients. A list
#' containing one entry, `scale`.
#' @param vertical_basis_function_variance_prior A prior object for the variance
#' of the vertical basis functions, given as a list with `shape` and
#' `rate` parameters. Defaults to an inverse gamma prior with 5% and 95%
#' percentiles of `1e-3 ^ 2` and `10 ^ 2` respectively
#' @param vertical_basis_function_boundary_knots The number of knots to place
#' before after after the domain in order to address boundary effects. Defaults
#' to three.
#'
#' @return A list containing the specifications for the GeoWarp variance model,
#' ready to be passed as an argument to \code{\link{geowarp_deviation_model}}
#' and \code{\link{geowarp_model}}.
#'
#' @examples
#' # Creating a GeoWarp variance model object with custom settings
#' variance_model <- geowarp_variance_model(
#'   fixed_formula = ~ depth,
#'   vertical_basis_function_delta = 0.5
#' )
#'
#' @seealso
#' \code{\link{geowarp_model}}, \code{\link{geowarp_deviation_model}},
#' \code{\link{geowarp_mean_model}}
#' @export
geowarp_variance_model <- function(
  fixed_formula = ~ 1,
  fixed_effect_mean = 0,
  fixed_effect_precision = 1 / 100,
  vertical_basis_functions = TRUE,
  vertical_basis_function_delta = 1,
  vertical_basis_function_length_scale_prior = list(scale = 1),
  vertical_basis_function_variance_prior = inverse_gamma_quantile_prior(
    1e-3 ^ 2,
    10 ^ 2
  ),
  vertical_basis_function_boundary_knots = 3L
) {
  list(
    fixed_formula = fixed_formula,
    fixed_effect_mean = fixed_effect_mean,
    fixed_effect_precision = fixed_effect_precision,
    vertical_basis_functions = vertical_basis_functions,
    vertical_basis_function_delta = vertical_basis_function_delta,
    vertical_basis_function_length_scale_prior = vertical_basis_function_length_scale_prior,
    vertical_basis_function_variance_prior = vertical_basis_function_variance_prior,
    vertical_basis_function_boundary_knots = vertical_basis_function_boundary_knots
  )
}

#' Create a GeoWarp Axial Warping Unit (AWU) Object
#'
#' Constructs a GeoWarp warping unit object for use in the deviation model
#' (see \code{\link{geowarp_deviation_model}}).
#'
#' Axial warping units operate on individual coordinates, and can be either
#' linear (\code{\link{geowarp_linear_awu}}) or non-linear using a Bernstein
#' polynomial (\code{\link{geowarp_bernstein_awu}}). A geometric warping unit
#' induces geometric anisotropy in the space, and is specified through
#' \code{\link{geowarp_geometric_warping_unit}}.
#'
#' @param scaling Numeric value specifying the scaling for the linear warping
#' unit. This is used to rescale the input coordinate, and the specified prior
#' applies after rescaling. This help with numerical issues by making the
#' inferred parameter closer to one.
#' @param prior A list giving the parameters for the gamma prior over the
#' axial warping unit parameters. This contains five entries: `type`, `shape`
#' `rate`, `lower`, and `upper`. The `type` parameter specifies the distribution
#' type, one of 'gamma', 'inv_uniform', or 'uniform'. The `shape` and `rate`
#' parameters specify the shape and rate parameters of the gamma distribution,
#' if chosen. The `lower` and `upper` parameters specify the lower and upper
#' bounds of the prior distribution; these can be set to `-Inf` or `Inf` as
#' needed.
#' @param order Numeric value specifying the order of the Bernstein polynomial.
#' Defaults to 20.
#' @param prior_shape Numeric value specifying the shape parameter for the
#' prior for the geometric warping unit. Defaults to 2. The value 1 implies a
#' uniform prior over the warping, while values greater than 1 imply shrinkage
#' towards no warping.
#'
#' @return A list representing the warping unit, ready to be passed as an
#' argument to \code{\link{geowarp_deviation_model}}.
#'
#' @examples
#' # Creating a GeoWarp linear AWU
#' linear_awu <- geowarp_linear_awu(scaling = 0.5)
#'
#' # Creating a GeoWarp Bernstein AWU
#' bernstein_awu <- geowarp_bernstein_awu(order = 10)
#'
#' # Creating a GeoWarp geometric warping unit
#' geometric_wu <- geowarp_geometric_warping_unit(prior_shape = 3)
#'
#' @seealso
#' \code{\link{geowarp_model}}, \code{\link{geowarp_deviation_model}},
#' \code{\link{geowarp_variance_model}}
#' @export
geowarp_linear_awu <- function(
  scaling = 1,
  prior = list(
    type = 'gamma',
    shape = 1.01,
    rate = 0.01,
    lower = 0,
    upper = Inf
  )
) {
  list(name = 'linear_awu', scaling = scaling, prior = prior)
}

#' @describeIn geowarp_linear_awu A warping unit that uses Bernstein
#' polynomials.
#' @export
geowarp_bernstein_awu <- function(
  order = 20,
  scaling = 1,
  prior = list(
    type = 'gamma',
    shape = 1.01,
    rate = 0.01,
    lower = 0,
    upper = Inf
  )
) {
  list(
    name = 'bernstein_awu',
    order = order,
    scaling = scaling,
    prior = prior
  )
}

#' @describeIn geowarp_linear_awu A geometric warping unit that specifies
#' geometric anisotropy.
#' @export
geowarp_geometric_warping_unit <- function(prior_shape = 2) {
  list(
    name = 'geometric_warping_unit',
    prior_shape = prior_shape
  )
}

#' Produce Matrix of Coordinates for a GeoWarp Model
#'
#' Extracts the coordinates specified in the given GeoWarp model from the input
#' data frame and returns them as a matrix.
#'
#' @param df A data frame containing the coordinates.
#' @param model A GeoWarp model object made using \code{\link{geowarp_model}}.
#'
#' @return A matrix containing the coordinates extracted from the data frame
#' according to the specifications in the GeoWarp model object. The columns
#' are ordered according to `model$horizontal_coordinates` followed by
#' `model$vertical_coordinate`.
#'
#' @examples
#' # Sample data frame
#' df <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6), z = c(7, 8, 9))
#'
#' # Sample GeoWarp model object
#' model <- geowarp_model(
#'   variable = 'value',
#'   horizontal_coordinates = c('x', 'y'),
#'   horizontal_domains = list(c(0, 10), c(0, 10)),
#'   vertical_coordinate = 'z',
#'   vertical_domain = c(0, 10)
#' )
#'
#' # Generate matrix of coordinates
#' coordinates_matrix <- model_coordinates(df, model)
#'
#' @seealso
#' \code{\link{geowarp_model}}, \code{\link{geowarp_deviation_model}},
#' \code{\link{geowarp_variance_model}}
#' @export
model_coordinates <- function(df, model) {
  output <- as.matrix(df[, c(model$horizontal_coordinates, model$vertical_coordinate)])
  dimnames(output) <- NULL
  output
}

.check_horizontal_warping <- function(warpings) {
  warping_names <- sapply(warpings, getElement, 'name')
  if (!all(head(warping_names, -1) == 'linear_awu')) {
    stop('non-linear warping currently not supported for horizontal dimensions')
  }
}
