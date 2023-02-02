#' @export
pcpt_optimise <- function(
  df,
  model = pcpt_model(),
  max_attempts = 20,
  best_of = 10,
  show_progress = FALSE,
  algorithm = 'BFGS',
  vecchia = 'auto',
  vecchia_n_parents = 20,
  vecchia_scaling = 'auto',
  vecchia_grouping_exponent = 2,
  parent_structure = get_parent_structure(
    df,
    model,
    vecchia_n_parents,
    scaling = vecchia_scaling
  ),
  ...
) {
  is_white <- model$deviation_model$name == 'white'
  if (vecchia == 'auto') {
    vecchia <- nrow(df) > 1000
  }
  vecchia <- is_white || vecchia

  log_debug('Converting inputs into format required for Stan')
  stan_data <- .to_stan_data(df, model)

  if (vecchia) {
    stan_data <- .augment_stan_data_with_vecchia(
      stan_data,
      vecchia_n_parents,
      parent_structure,
      model,
      vecchia_grouping_exponent
    )
  }

  log_debug('Running optimization')
  model_name <- if (is_white) 'white' else sprintf(
    '%s_%s',
    model$deviation_model$name,
    if (vecchia) 'vecchia' else 'exact'
  )
  stan_fit <- .optimizing_best_of(
    max_attempts = max_attempts,
    best_of = best_of,
    show_progress = show_progress,
    object = stanmodels[[model_name]],
    data = stan_data,
    as_vector = FALSE,
    hessian = FALSE,
    algorithm = algorithm,
    ...
  )
  parameters <- stan_fit$par
  if (!is_white) {
    parameters$alpha_beta <- parameters$alpha_beta_hat
    parameters$alpha_beta_hat <- NULL
    parameters$log_marginal <- NULL
  }
  structure(
    list(
      model = model,
      input_df = df,
      parameters = parameters,
      stan_fit = stan_fit
    ),
    class = 'pcpt_fit'
  )
}

#' @export
pcpt_sample <- function(
  df,
  model = pcpt_model(),
  vecchia = 'auto',
  vecchia_n_parents = 20,
  vecchia_scaling = 'auto',
  vecchia_grouping_exponent = 2,
  parent_structure = get_parent_structure(
    df,
    model,
    vecchia_n_parents,
    scaling = vecchia_scaling
  ),
  control = list(
    metric = 'dense_e'
  ),
  ...
) {
  is_white <- model$deviation_model$name == 'white'
  if (vecchia == 'auto') {
    vecchia <- nrow(df) > 1000
  }
  vecchia <- is_white || vecchia

  log_debug('Converting inputs into format required for Stan')
  stan_data <- .to_stan_data(df, model)

  if (vecchia) {
    stan_data <- .augment_stan_data_with_vecchia(
      stan_data,
      vecchia_n_parents,
      parent_structure,
      model,
      vecchia_grouping_exponent
    )
  }

  log_debug('Running optimization')
  model_name <- if (is_white) 'white' else sprintf(
    '%s_%s',
    model$deviation_model$name,
    if (vecchia) 'vecchia' else 'exact'
  )
  stan_fit <- rstan::sampling(
    object = stanmodels[[model_name]],
    data = stan_data,
    ...
  )

  structure(
    list(
      model = model,
      input_df = df,
      stan_fit = stan_fit
    ),
    class = 'pcpt_fit'
  )
}

# #' @export
# print.pcpt_fit <- function(x, ...) {
#   cat('- Model:\n')
#   print(x$model)
#   cat('- Estimated parameters:\n')
#   if (is(x$stan_fit, 'stanfit')) {
#     print(summary(x$stan_fit))
#   } else {
#     print(x$stan_fit$par)
#   }
#   invisible(x)
# }

.to_stan_data <- function(
  df,
  model
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
      output$scaling <- c(rep(1, output$D - 1), model$deviation_model$warping$scaling)
    } else {
      output$scaling <- sapply(model$deviation_model$warpings, getElement, 'scaling')
    }
  }

  output$y <- df[[model$variable]]

  if (is_vertical_only) {
    ordering <- do.call(order, lapply(
      c(model$horizontal_coordinates, model$vertical_coordinate),
      function(name) df[[name]]
    ))
    stopifnot(!is.unsorted(ordering))

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
  output$tau_squared_mean_random_scale <- model$mean_model$vertical_basis_function_variance_prior$scale

  # Deviation model
  if (!is_white) {
    vertical_warping <- if (is_vertical_only) {
      model$deviation_model$warping
    } else {
      model$deviation_model$warpings[[output$D]]
    }
    if (vertical_warping$name == 'linear') {
      output$X_deviation_warping <- cbind(
        vertical_warping$scaling * output$x[, output$D]
      )
      output$P_deviation_warping <- 1L
    } else {
      output$X_deviation_warping <- bernstein_warping_design_matrix(
        vertical_warping$scaling * output$x[, output$D],
        vertical_warping$order,
        vertical_warping$domain
      )
      output$P_deviation_warping <- ncol(output$X_deviation_warping)
    }
  }

  deviation_parts <- .random_effects_design_matrices(
    df,
    model$vertical_coordinate,
    model$deviation_model$variance_model
  )
  output$X_deviation_fixed <- deviation_parts$X_fixed
  output$X_deviation_random <- deviation_parts$X_random
  output$P_deviation_fixed <- ncol(output$X_deviation_fixed)
  output$P_deviation_random <- ncol(output$X_deviation_random)
  output$eta_deviation_mean <- array(vctrs::vec_recycle(
    model$deviation_model$variance_model$fixed_effect_mean,
    output$P_deviation_fixed
  ), dim = output$P_deviation_fixed)
  if (is.matrix(model$deviation_model$variance_model$fixed_effect_precision)) {
    output$eta_deviation_precision <- model$deviation_model$variance_model$fixed_effect_precision
  } else {
    output$eta_deviation_precision <- diag(
      x = vctrs::vec_recycle(
        model$deviation_model$variance_model$fixed_effect_precision,
        output$P_deviation_fixed
      ),
      nrow = output$P_deviation_fixed,
      ncol = output$P_deviation_fixed
    )
  }
  output$delta_deviation_random <- model$deviation_model$variance_model$vertical_basis_function_delta
  output$ell_deviation_random_scale <- model$deviation_model$variance_model$vertical_basis_function_length_scale_prior$scale
  output$tau_squared_deviation_random_scale <- model$deviation_model$variance_model$vertical_basis_function_variance_prior$scale

  if (!is_white) {
    if (is_vertical_only) {
      output$gamma_deviation_a <- model$deviation_model$warping$prior$shape
      output$gamma_deviation_b <- model$deviation_model$warping$prior$rate
    } else {
      warping_priors <- lapply(model$deviation_model$warpings, getElement, 'prior')
      output$gamma_deviation_a <- sapply(warping_priors, getElement, 'shape')
      output$gamma_deviation_b <- sapply(warping_priors, getElement, 'rate')
    }
  }

  if (model$deviation_model$name == 'anisotropic') {
    output$L_deviation_shape <- model$deviation_model$anisotropy_shape
  }

  output$sigma_squared_nugget_a <- model$nugget_prior$shape
  output$sigma_squared_nugget_b <- model$nugget_prior$rate

  output
}

.augment_stan_data_with_vecchia <- function(stan_data, n_parents, parent_structure, model, vecchia_grouping_exponent) {
  if (model$deviation_model$name == 'vertical_only') {
    stan_data$N_blocks <- 0L
    stan_data$block_indices <- integer(0)
    stan_data$block_last_index <- integer(0)
    stan_data$block_N_responses <- integer(0)
    for (i in 1 : stan_data$N_individuals) {
      if (i == stan_data$N_individuals) {
        end <- stan_data$N
      } else {
        end <- stan_data$start[i + 1] - 1
      }

      N_i <- end - stan_data$start[i] + 1
      N_blocks_i <- ceiling(N_i / n_parents)

      for (j in 1 : N_blocks_i) {
        stan_data$N_blocks <- stan_data$N_blocks + 1
        response_start <- (j - 1) * n_parents + 1
        response_end <- min(N_i, j * n_parents)
        parent_start <- max(1, response_start - n_parents)
        parent_end <- max(1, response_start - 1)

        stan_data$block_indices <- c(
          stan_data$block_indices,
          stan_data$start[i] + (parent_start : response_end) - 1L
        )
        stan_data$block_last_index <- c(stan_data$block_last_index, length(stan_data$block_indices))
        stan_data$block_N_responses <- c(stan_data$block_N_responses, response_end - response_start + 1L)
      }
    }
    stan_data$N_indices <- length(stan_data$block_indices)
  } else if (model$deviation_model$name == 'white') {
    stan_data$block_indices <- seq_len(stan_data$N)
    stan_data$block_last_index <- tail(seq(0, stan_data$N, by = n_parents), -1)
    if (tail(stan_data$block_last_index, 1) != stan_data$N) {
      stan_data$block_last_index <- c(
        stan_data$block_last_index,
        stan_data$N
      )
    }
    stan_data$block_N_responses <- diff(c(0, stan_data$block_last_index))
    stan_data$N_indices <- length(stan_data$block_indices)
    stan_data$N_blocks <- length(stan_data$block_N_responses)
  } else {
    stan_data$y <- stan_data$y[parent_structure$ordering]
    for (name in c('x', 'X_mean_fixed', 'X_mean_random', 'X_deviation_warping', 'X_deviation_fixed', 'X_deviation_random')) {
      stan_data[[name]] <- stan_data[[name]][parent_structure$ordering, , drop = FALSE]
    }
    vecchia_groups <- GpGp::group_obs(parent_structure$parents, exponent = vecchia_grouping_exponent)
    current_block_index <- 1L
    current_response_index <- 1L
    for (i in seq_along(vecchia_groups$last_ind_of_block)) {
      block_i <- current_block_index : vecchia_groups$last_ind_of_block[i]
      block_indices <- vecchia_groups$all_inds[block_i]
      response_indices <- vecchia_groups$local_resp_inds[current_response_index : vecchia_groups$last_resp_of_block[i]]
      vecchia_groups$all_inds[block_i] <- c(
        block_indices[-response_indices],
        block_indices[response_indices]
      )
      current_block_index <- vecchia_groups$last_ind_of_block[i] + 1L
      current_response_index <- vecchia_groups$last_resp_of_block[i] + 1L
    }
    stan_data$N_indices <- length(vecchia_groups$all_inds)
    stan_data$N_blocks <- length(vecchia_groups$last_ind_of_block)
    stan_data$block_indices <- vecchia_groups$all_inds
    stan_data$block_last_index <- vecchia_groups$last_ind_of_block
    stan_data$block_N_responses <- c(vecchia_groups$last_resp_of_block[1], diff(vecchia_groups$last_resp_of_block))
  }
  stan_data
}

.random_effects_design_matrices <- function(df, vertical_coordinate, model) {
  output <- NULL
  output$X_fixed <- model.matrix(model$fixed_formula, df)
  if (model$vertical_basis_functions) {
    x <- df[[vertical_coordinate]]

    shift_value <- function(z) {
      (z - model$vertical_basis_function_domain[1]) / model$vertical_basis_function_delta
    }

    centres <- seq(
      model$vertical_basis_function_delta * (floor(shift_value(
        model$vertical_basis_function_domain[1]
      )) - 3),
      model$vertical_basis_function_delta * (ceiling(shift_value(
        model$vertical_basis_function_domain[2]
      )) + 3),
      by = model$vertical_basis_function_delta
    )
    output$X_random <- splines::splineDesign(
      knots = centres,
      x = x
    )
  } else {
    output$X_random <- matrix(0, nrow = nrow(df), ncol = 0)
  }
  output
}

.optimizing_best_of <- function(max_attempts, best_of, show_progress = FALSE, ...) {
  good_attempts <- 0
  best_result <- list(value = -Inf)
  attempts <- list()
  for (attempt in seq_len(max_attempts)) {
    if (show_progress) {
      cat('\rrunning', attempt, '/', max_attempts, '| completed', good_attempts, '/', best_of, '| best', best_result$value)
    }
    result <- tryCatch({
      rstan::optimizing(...)
    }, error = function(e) {
      print(e)
      NULL
    })
    if (!is.null(result)) {
      attempts <- c(attempts, list(result))
    }
    if (is.null(result)) next

    good_attempts <- good_attempts + 1
    if (result$value > best_result$value) {
      best_result <- result
    }
    if (good_attempts == best_of) break
  }
  cat('\r')
  if (attempt == max_attempts) stop('Max attempts exceeded')

  best_result$attempts <- attempts
  best_result
}
