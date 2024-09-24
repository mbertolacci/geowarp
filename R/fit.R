#' Estimate the GeoWarp Model Parameters Through Optimisation
#'
#' This function estimates the parameters of a specified GeoWarp model by using
#' an optimisation algorithm to find the posterior mode.
#'
#' @param df Data frame containing the observations.
#' @param model A GeoWarp model object specifying the model structure and
#' priors,created using \code{\link{geowarp_model}}.
#' @param max_attempts Maximum number of optimisation attempts. Defaults to 20.
#' @param best_of Number of times to run the optimisation algorithm. The
#' algorithm is run repeatedly with different starting values in order to help
#' find the global posterior mode.
#' @param method Optimisation algorithm to use, one of either 'stan_trust_optim'
#' or 'optimizing'. The first, the default, uses the
#' \code{\link{stan_trust_optim}} function, while the second uses
#' \code{\link[rstan]{optimizing}}.
#' @param vecchia Whether to use the Vecchia approximation. Can be 'auto',
#' `TRUE`, or `FALSE`. If 'auto', the Vecchia approximation is used if there
#' are more than 1000 observations.
#' @param n_parents Number of parent locations for the Vecchia approximation.
#' Defaults to 20.
#' @param parent_structure The parent structure used for the Vecchia
#' approximation; see \code{\link{vecchia_parent_structure}}.
#' @param grouping Grouping structure for the Vecchia approximation; see
#' \code{\link{vecchia_grouping}}. The default should be fine for most cases.
#' @param threads Number of threads to use for parallelisation, taken from
#' `getOption('geowarp.threads')` by default, which is itself set to 1 by
#' default. The special value -1 picks a number of threads based on the number
#' of cores in the system).
#' @param grain_size Used to control the parallelisation. The default, 100,
#' should be okay. On systems with few cores, a larger value might be better.
#' @param ... Additional parameters to be passed to the chosen optimisation
#' routine in \code{method}.
#'
#' @return A 'geowarp_fit' object that estimated parameters, suitable for use
#' with \code{\link{predict.geowarp_fit}} and other methods.
#'
#' @examples
#' # Using the geowarp_optimise function
#' fit <- geowarp_optimise(
#'   df = my_data,
#'   model = my_model
#' )
#'
#' @export
geowarp_optimise <- function(
  df,
  model = geowarp_model(),
  max_attempts = 20,
  best_of = 5,
  method = c('stan_nlminb', 'stan_nlm', 'stan_trust_optim', 'optimizing'),
  vecchia = 'auto',
  n_parents = 20,
  parent_structure = vecchia_parent_structure(
    df,
    model,
    n_parents,
  ),
  grouping = vecchia_grouping(parent_structure),
  threads = getOption('geowarp.threads'),
  grain_size = 100L,
  ...
) {
  method <- match.arg(method)

  .local_tbb_threads(threads)
  .validate_data(df, model)

  is_white <- model$deviation_model$name == 'white'
  if (vecchia == 'auto') {
    vecchia <- nrow(df) > 1000
  }
  if (is_white) {
    vecchia <- FALSE
  }

  log_debug('Converting inputs into format required for Stan')
  stan_data <- geowarp_stan_data(
    df,
    model,
    threads,
    grain_size,
    if (vecchia) parent_structure else NULL,
    if (vecchia) grouping else NULL
  )

  log_debug('Running optimisation')
  stan_fit <- if (method != 'optimizing') {
    stan_optimizing_best_of(
      object = stanmodels[[model$deviation_model$name]],
      data = stan_data,
      .max_attempts = max_attempts,
      .best_of = best_of,
      .method = match.fun(method),
      ...
    )
  } else {
    stan_optimizing_best_of(
      object = stanmodels[[model$deviation_model$name]],
      data = stan_data,
      as_vector = FALSE,
      hessian = FALSE,
      .max_attempts = max_attempts,
      .best_of = best_of,
      .method = rstan::optimizing,
      ...
    )
  }
  parameters <- stan_fit$par
  parameters$alpha_beta <- parameters$alpha_beta_hat
  parameters$alpha_beta_hat <- NULL
  parameters$log_marginal <- NULL
  output <- structure(
    list(
      model = model,
      observed_df = df,
      parameters = parameters,
      stan_fit = stan_fit
    ),
    class = 'geowarp_fit'
  )
  if (vecchia) {
    output$parent_structure <- parent_structure
    output$grouping <- grouping
  }
  output
}

#' GeoWarp Stan Model
#'
#' These functions allow access to the underlying Stan implementation of
#' GeoWarp. They allow access to the \link[rstan]{stan_model} object used to fit
#' GeoWarp, and can prepare the input data required to use the object.
#'
#' @param name Name of the model to use. Can be 'full', 'vertical_only', or
#' 'white'. The full model is the default.
#' @param df The input data to prepare for Stan.
#' @param model The GeoWarp model object, created using
#' \code{\link{geowarp_model}}.
#' @param parent_structure The parent structure used for the Vecchia
#' approximation; see \code{\link{vecchia_parent_structure}}.
#' @param grouping Grouping structure for the Vecchia approximation; see
#' \code{\link{vecchia_grouping}}. The default should be fine for most cases.
#' @param white_block_size Block size to use for the white noise model. This
#' affects performance but not the results; the default should be fine.
#'
#' @return A \link[rstan]{stan_model} object.
#'
#' @seealso \code{\link{geowarp_optimise}}
#'
#' @export
geowarp_stan_model <- function(name = c('full', 'vertical_only', 'white')) {
  name <- match.arg(name)
  stanmodels[[name]]
}

#' @describeIn geowarp_stan_model Prepare data for Stan
#' @export
geowarp_stan_data <- function(
  df,
  model,
  threads = 1L,
  grain_size = 100L,
  parent_structure,
  grouping,
  white_block_size = min(nrow(df) - 1L, 30)
) {
  is_white <- model$deviation_model$name == 'white'
  is_vertical_only <- model$deviation_model$name == 'vertical_only'

  output <- list()

  # Data
  output$N <- nrow(df)
  output$D <- length(c(model$horizontal_coordinates, model$vertical_coordinate))
  output$x <- as.matrix(
    df[, c(model$horizontal_coordinates, model$vertical_coordinate)]
  )
  if (!is_white) {
    if (is_vertical_only) {
      output$scaling <- c(
        rep(1, output$D - 1),
        model$deviation_model$axial_warping_unit$scaling
      )
    } else {
      output$scaling <- sapply(
        model$deviation_model$axial_warping_units,
        getElement, 'scaling'
      )[model$deviation_model$axial_warping_unit_mapping]
    }
  }

  output$y <- df[[model$variable]]

  if (is_vertical_only) {
    coordinates <- df %>%
      group_by(across(model$horizontal_coordinates)) %>%
      summarise(n = n())
    output$N_individuals <- nrow(coordinates)
    output$start <- c(1, 1 + head(cumsum(coordinates$n), -1))
  }

  # Mean model
  mean_parts <- .random_effects_design_matrices(
    df,
    model$vertical_coordinate,
    model$vertical_domain,
    model$mean_model
  )
  output$X_mean_fixed <- mean_parts$X_fixed
  output$X_mean_random <- mean_parts$X_random
  output$P_mean_fixed <- ncol(mean_parts$X_fixed)
  output$P_mean_random <- ncol(mean_parts$X_random)
  output$alpha_mean <- array(vctrs::vec_recycle(
    model$mean_model$fixed_effect_mean,
    output$P_mean_fixed
  ), dim = output$P_mean_fixed)
  if (is.matrix(model$mean_model$fixed_effect_precision)) {
    output$alpha_precision <- model$mean_model$fixed_effect_precision
  } else {
    output$alpha_precision <- diag(
      x = vctrs::vec_recycle(
        model$mean_model$fixed_effect_precision,
        output$P_mean_fixed
      ),
      nrow = output$P_mean_fixed,
      ncol = output$P_mean_fixed
    )
  }
  output$tau_squared_mean_random_a <- model$mean_model$vertical_basis_function_variance_prior$shape
  output$tau_squared_mean_random_b <- model$mean_model$vertical_basis_function_variance_prior$rate

  # Deviation model
  if (!is_white) {
    output$smoothness <- switch(
      model$deviation_model$covariance_function,
      exponential = 0.5,
      matern15 = 1.5,
      stop('Covariance function not supported')
    )
    vertical_warping <- if (is_vertical_only) {
      model$deviation_model$axial_warping_unit
    } else {
      tail(model$deviation_model$axial_warping_units, 1)[[1]]
    }
    if (vertical_warping$name == 'linear_awu') {
      output$X_deviation_warping <- cbind(output$x[, output$D])
      output$P_deviation_warping <- 1L
    } else {
      output$X_deviation_warping <- bernstein_warping_design_matrix(
        output$x[, output$D],
        vertical_warping$order,
        model$vertical_domain
      )
      output$P_deviation_warping <- ncol(output$X_deviation_warping)
    }
    output$X_deviation_warping <- vertical_warping$scaling * output$X_deviation_warping
  }

  variance_model <- model$deviation_model$variance_model
  deviation_parts <- .random_effects_design_matrices(
    df,
    model$vertical_coordinate,
    model$vertical_domain,
    variance_model
  )
  output$X_deviation_fixed <- deviation_parts$X_fixed
  output$X_deviation_random <- deviation_parts$X_random
  output$P_deviation_fixed <- ncol(output$X_deviation_fixed)
  output$P_deviation_random <- ncol(output$X_deviation_random)
  output$eta_deviation_mean <- array(vctrs::vec_recycle(
    variance_model$fixed_effect_mean,
    output$P_deviation_fixed
  ), dim = output$P_deviation_fixed)
  if (is.matrix(variance_model$fixed_effect_precision)) {
    output$eta_deviation_precision <- variance_model$fixed_effect_precision
  } else {
    output$eta_deviation_precision <- diag(
      x = vctrs::vec_recycle(
        variance_model$fixed_effect_precision,
        output$P_deviation_fixed
      ),
      nrow = output$P_deviation_fixed,
      ncol = output$P_deviation_fixed
    )
  }
  output$delta_deviation_random <- variance_model$vertical_basis_function_delta
  output$ell_deviation_random_scale <- variance_model$vertical_basis_function_length_scale_prior$scale
  output$tau_squared_deviation_random_a <- variance_model$vertical_basis_function_variance_prior$shape
  output$tau_squared_deviation_random_b <- variance_model$vertical_basis_function_variance_prior$rate

  if (!is_white) {
    if (is_vertical_only) {
      output$gamma_deviation_prior_type <- model$deviation_model$axial_warping_unit$prior$type
      output$gamma_deviation_a <- model$deviation_model$axial_warping_unit$prior$shape
      output$gamma_deviation_b <- model$deviation_model$axial_warping_unit$prior$rate
      output$gamma_deviation_lower <- model$deviation_model$axial_warping_unit$prior$lower
      output$gamma_deviation_upper <- model$deviation_model$axial_warping_unit$prior$upper
    } else {
      warping_priors <- lapply(model$deviation_model$axial_warping_units, getElement, 'prior')
      output$gamma_deviation_prior_type <- sapply(warping_priors, getElement, 'type')
      output$gamma_deviation_a <- sapply(warping_priors, getElement, 'shape')
      output$gamma_deviation_b <- sapply(warping_priors, getElement, 'rate')
      output$gamma_deviation_lower <- sapply(warping_priors, getElement, 'lower')
      output$gamma_deviation_upper <- sapply(warping_priors, getElement, 'upper')
      output$D_horizontal_warpings <- length(warping_priors) - 1L
      output$axial_warping_unit_mapping <- model$deviation_model$axial_warping_unit_mapping
    }

    output$gamma_deviation_prior_type <- c(
      'gamma' = 1,
      'inv_uniform' = 2,
      'uniform' = 3
    )[output$gamma_deviation_prior_type]
  }

  if (!is.null(model$deviation_model$geometric_warping_unit)) {
    output$L_deviation_shape <- model$deviation_model$geometric_warping_unit$prior_shape
    output$D_geometric <- output$D
  } else {
    output$L_deviation_shape <- 0
    output$D_geometric <- 0L
  }

  output$sigma_squared_nugget_a <- model$nugget_prior$shape
  output$sigma_squared_nugget_b <- model$nugget_prior$rate

  ## Add Vecchia bits as needed
  has_grouping <- !missing(grouping) && !is.null(grouping)
  if (!has_grouping && !is_white) {
    # Create one block with all indices in it
    parent_structure <- list(
      observed_ordering = seq_len(output$N)
    )
    grouping <- list(
      N_indices = output$N,
      N_blocks = 1L,
      block_indices = seq_len(output$N),
      block_last_index = output$N,
      block_N_responses = output$N
    )
  }

  if (is_white) {
    output$block_indices <- seq_len(output$N)
    output$block_last_index <- tail(seq(0, output$N, by = white_block_size), -1)
    if (tail(output$block_last_index, 1) != output$N) {
      output$block_last_index <- c(
        output$block_last_index,
        output$N
      )
    }
    output$block_N_responses <- diff(c(0, output$block_last_index))
    output$N_indices <- length(output$block_indices)
    output$N_blocks <- length(output$block_N_responses)
  } else {
    output$y <- output$y[parent_structure$observed_ordering]
    for (name in c('x', 'X_mean_fixed', 'X_mean_random', 'X_deviation_warping', 'X_deviation_fixed', 'X_deviation_random')) {
      output[[name]] <- output[[name]][parent_structure$observed_ordering, , drop = FALSE]
    }
    output$N_indices <- grouping$N_indices
    output$N_blocks <- grouping$N_blocks
    output$block_indices <- array(as.integer(grouping$block_indices))
    output$block_last_index <- array(as.integer(grouping$block_last_index))
    output$block_N_responses <- array(as.integer(grouping$block_N_responses))
  }

  if (!is_white) {
    output$use_parallel <- threads == -1L || threads > 1L
    output$grain_size <- grain_size
  }

  output
}

.random_effects_design_matrices <- function(df, vertical_coordinate, vertical_domain, model) {
  output <- NULL
  output$X_fixed <- model.matrix(model$fixed_formula, df)
  if (model$vertical_basis_functions) {
    x <- df[[vertical_coordinate]]

    shift_value <- function(z) {
      (z - vertical_domain[1]) / model$vertical_basis_function_delta
    }

    centres <- seq(
      model$vertical_basis_function_delta * (floor(shift_value(
        vertical_domain[1]
      )) - model$vertical_basis_function_boundary_knots),
      model$vertical_basis_function_delta * (ceiling(shift_value(
        vertical_domain[2]
      )) + model$vertical_basis_function_boundary_knots),
      by = model$vertical_basis_function_delta
    )
    output$X_random <- splines::splineDesign(
      knots = centres,
      x = x,
      outer.ok = TRUE
    )
  } else {
    output$X_random <- matrix(0, nrow = nrow(df), ncol = 0)
  }
  output
}

.validate_data <- function(df, model) {
  if (
    anyDuplicated(rbind(as.matrix(
      df[, c(model$horizontal_coordinates, model$vertical_coordinate)]
    )))
  ) {
    stop('Duplicate coordinates are not supported')
  }
  coordinate_names <- c(model$horizontal_coordinates, model$vertical_coordinate)
  coordinate_domains <- c(
    model$horizontal_domains,
    list(model$vertical_domain)
  )
  for (i in seq_along(coordinate_names)) {
    coordinate_range <- range(df[[coordinate_names[i]]])
    if (
      coordinate_range[1] < coordinate_domains[[i]][1]
      || coordinate_range[2] > coordinate_domains[[i]][2]
    ) {
      stop(sprintf(
        'Coordinate %s is outside of the domain',
        coordinate_names[i]
      ))
    }
  }
}
