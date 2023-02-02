#' @export
mean_profile <- function(
  object,
  output_df = object$input_df
) {
  stan_data <- .to_stan_data(output_df, object$model)

  as.vector(cbind(
    stan_data$X_mean_fixed,
    stan_data$X_mean_random
  ) %*% object$parameters$alpha_beta)
}

#' @export
predict.pcpt_fit <- function(
  object,
  output_df = object$input_df,
  input_df = object$input_df,
  include_mean = TRUE,
  include_precision_V = FALSE,
  include_samples = FALSE,
  n_samples = 500,
  vecchia = 'auto',
  vecchia_n_parents = 50,
  parent_structure = get_parent_structure(
    list(input_df, output_df),
    object$model,
    vecchia_n_parents,
    scaling = if (object$model$deviation_model$name == 'vertical_only') {
      c(rep(1, length(object$model$horizontal_coordinates)), 1e-6)
    } else {
      'auto'
    },
    ordering = if (object$model$deviation_model$name == 'vertical_only') {
      function(x) {
        coordinates <- lapply(seq_len(ncol(x)), function(i) {
          x[, i]
        })
        do.call(order, coordinates)
      }
    } else {
      GpGp::order_maxmin
    }
  ),
  ...
) {
  if (vecchia == 'auto') {
    vecchia <- (nrow(input_df) + nrow(output_df)) > 1000
  }

  stan_data <- lapply(
    list(input_df, output_df),
    .to_stan_data,
    object$model
  )

  if (object$model$deviation_model$name == 'white') {
    prediction_distribution <- .prediction_distribution_white(
      stan_data[[1]],
      stan_data[[2]],
      object
    )
  } else if (!vecchia) {
    prediction_distribution <- .prediction_distribution_exact(
      stan_data[[1]],
      stan_data[[2]],
      object
    )
  } else {
    prediction_distribution <- .prediction_distribution_vecchia(
      stan_data[[1]],
      stan_data[[2]],
      object,
      parent_structure
    )
  }
  output <- list()
  if (include_mean) {
    output$mean <- prediction_distribution$mean
  }
  if (include_precision_V) {
    output$precision_V <- prediction_distribution$precision_V
    output$ordering <- prediction_distribution$ordering
  }
  if (include_samples) {
    z <- matrix(rnorm(nrow(output_df) * n_samples), nrow = nrow(output_df))
    if (!is(prediction_distribution$precision_V, 'sparseMatrix')) {
      y_tilde <- backsolve(
        prediction_distribution$precision_V,
        z,
        transpose = TRUE
      )
    } else {
      y_tilde <- as.matrix(solve(
        t(prediction_distribution$precision_V),
        z
      ))
    }
    output$samples <- prediction_distribution$mean + y_tilde[prediction_distribution$ordering, ]
  }
  output
}

.prediction_distribution_exact <- function(
  input_stan_data,
  output_stan_data,
  fit
) {
  rbind_matrices <- function(name) {
    rbind(input_stan_data[[name]], output_stan_data[[name]])
  }
  Sigma_all <- .covariance_matrix_internal(
    rbind_matrices('x'),
    rbind_matrices('X_deviation_fixed'),
    rbind_matrices('X_deviation_random'),
    fit
  )

  n_input <- nrow(input_stan_data$x)
  n_output <- nrow(output_stan_data$x)
  input_indices <- seq_len(n_input)
  output_indices <- n_input + seq_len(n_output)

  Sigma_input <- Sigma_all[input_indices, input_indices]
  Sigma_output <- Sigma_all[output_indices, output_indices]
  Sigma_cross <- Sigma_all[input_indices, output_indices]

  mu_input <- as.vector(cbind(
    input_stan_data$X_mean_fixed,
    input_stan_data$X_mean_random
  ) %*% fit$parameters$alpha_beta)
  mu_output <- as.vector(cbind(
    output_stan_data$X_mean_fixed,
    output_stan_data$X_mean_random
  ) %*% fit$parameters$alpha_beta)

  chol_Sigma_input <- chol(Sigma_input)
  mu_prediction <- mu_output + as.vector(crossprod(Sigma_cross, .chol_solve(
    chol_Sigma_input,
    input_stan_data$y - mu_input
  )))
  Sigma_prediction <- (
    Sigma_output
    - crossprod(Sigma_cross, .chol_solve(chol_Sigma_input, Sigma_cross))
  )
  Q_prediction <- chol2inv(chol(Sigma_prediction))
  rev_matrix <- function(x) {
    x[rev(seq_len(nrow(x))), rev(seq_len(ncol(x)))]
  }
  list(
    mean = mu_prediction,
    precision_V = rev_matrix(t(chol(rev_matrix(Q_prediction)))),
    ordering = seq_len(n_output)
  )
}

.prediction_distribution_vecchia <- function(
  input_stan_data,
  output_stan_data,
  fit,
  parent_structure
) {
  rbind_matrices <- function(name) {
    rbind(
      input_stan_data[[name]][parent_structure$ordering[[1]], , drop = FALSE],
      output_stan_data[[name]][parent_structure$ordering[[2]], , drop = FALSE]
    )
  }
  x_all <- rbind_matrices('x')
  X_mean_all <- cbind(
    rbind_matrices('X_mean_fixed'),
    rbind_matrices('X_mean_random')
  )
  X_deviation_fixed_all <- rbind_matrices('X_deviation_fixed')
  X_deviation_random_all <- rbind_matrices('X_deviation_random')

  n_input <- nrow(input_stan_data$x)
  n_output <- nrow(output_stan_data$x)
  input_indices <- seq_len(n_input)
  output_indices <- n_input + seq_len(n_output)

  U <- .create_vecchia_U(
    x_all,
    X_deviation_fixed_all,
    X_deviation_random_all,
    fit,
    parent_structure
  )
  V <- U[output_indices, output_indices]

  mu_hat <- as.vector(X_mean_all %*% fit$parameters$alpha_beta)
  prediction_mean <- (
    mu_hat[output_indices]
    - as.vector(solve(t(V), solve(V,
      U[output_indices, ] %*%
      crossprod(
        U[input_indices, ],
        input_stan_data$y[parent_structure$ordering[[1]]]
        - mu_hat[input_indices]
      )
    )))
  )
  reverse_ordering <- invPerm(parent_structure$ordering[[2]])
  list(
    mean = prediction_mean[reverse_ordering],
    precision_V = V,
    ordering = reverse_ordering
  )
}

.prediction_distribution_white <- function(
  input_stan_data,
  output_stan_data,
  fit
) {
  prediction_mean <- as.vector(cbind(
    output_stan_data$X_mean_fixed,
    output_stan_data$X_mean_random
  ) %*% fit$parameters$alpha_beta)
  prediction_sd <- exp(as.vector(0.5 * (
    output_stan_data$X_deviation_fixed %*% fit$parameters$eta_deviation
    + output_stan_data$X_deviation_random %*% fit$parameters$zeta_deviation
  )))

  precision_V <- Diagonal(x = 1 / prediction_sd)
  list(
    mean = prediction_mean,
    precision_V = precision_V,
    ordering = seq_len(output_stan_data$N)
  )
}
