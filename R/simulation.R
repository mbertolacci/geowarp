#' @export
simulate_pcpts <- function(
  output_df,
  model,
  parameters,
  n_samples = 1,
  vecchia = 'auto',
  vecchia_n_parents = 50,
  parent_structure = vecchia_parent_structure(
    output_df,
    model,
    vecchia_n_parents
  )
) {
  if (vecchia == 'auto') {
    vecchia <- nrow(output_df) > 1000
  }
  stan_data <- .to_stan_data(output_df, model)
  output_df[[model$variable]] <- if (!vecchia) {
    .simulate_pcpts_exact(
      stan_data,
      model,
      parameters,
      n_samples
    )
  } else {
    .simulate_pcpts_vecchia(
      stan_data,
      model,
      parameters,
      parent_structure,
      n_samples
    )
  }
  output_df
}

.simulate_pcpts_exact <- function(
  stan_data,
  model,
  parameters,
  n_samples
) {
  Sigma <- .covariance_matrix_internal(
    stan_data$x,
    stan_data$X_deviation_fixed,
    stan_data$X_deviation_random,
    model = model,
    parameters = parameters
  )
  R <- chol(Sigma)

  output <- as.vector(cbind(
    stan_data$X_mean_fixed,
    stan_data$X_mean_random
  ) %*% parameters$alpha_beta) + crossprod(
    R,
    matrix(rnorm(n_samples * nrow(Sigma)), ncol = n_samples)
  )
  if (n_samples == 1) as.vector(output)
  else output
}

.simulate_pcpts_vecchia <- function(
  stan_data,
  model,
  parameters,
  parent_structure,
  n_samples
) {
  U <- .create_vecchia_U(
    stan_data$x[parent_structure$ordering, , drop = FALSE],
    stan_data$X_deviation_fixed[parent_structure$ordering, , drop = FALSE],
    stan_data$X_deviation_random[parent_structure$ordering, , drop = FALSE],
    parent_structure = parent_structure,
    model = model,
    parameters = parameters,
  )

  output <- as.vector(cbind(
    stan_data$X_mean_fixed,
    stan_data$X_mean_random
  )[parent_structure$ordering, , drop = FALSE] %*% parameters$alpha_beta) + solve(
    t(U),
    matrix(rnorm(n_samples * nrow(U)), ncol = n_samples)
  )
  output[Matrix::invPerm(parent_structure$ordering), ]
}
