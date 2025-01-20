#' Simulate from a GeoWarp Fit
#'
#' This function performs unconditional simulation from a fitted GeoWarp model.
#' The simulation conditions only on the estimated parameters and not on the
#' data. To simulate conditioned on data, use \code{\link{predict.geowarp_fit}}.
#'
#' @param fit GeoWarp fit object made using \code{\link{geowarp_optimise}}. This
#' is optional if `df`, `model`, and `parameters` are given directly.
#' @param df A data frame containing the locations to simulated (defaults to the
#' observation locations)
#' @param model GeoWarp model object used for simulation (default uses the model
#' from the `fit` object).
#' @param parameters Parameters used for the simulation (default uses the
#' parameters from the `fit` object).
#' @param n_samples The number of samples to simulate (default is 1).
#' @param nugget Logical indicating whether to include a nugget effect (default
#' is TRUE).
#' @param vecchia Either 'auto', TRUE, or FALSE. If 'auto', vecchia
#' approximation is used for n > 1000 (default is 'auto').
#' @param vecchia_n_parents Number of parents in the Vecchia approximation
#' (default is 50).
#' @param parent_structure The parent structure for Vecchia approximation
#' constructed using \code{\link{vecchia_parent_structure}}. The default parent
#' structure assumes that the data are vertically dense, which is not
#' appropriate if they are uniform or on a grid; see the documentation for that
#' function.
#' @param threads Number of threads to use for parallelisation, taken from
#' `getOption('geowarp.threads')` by default, which is itself set to 1 by
#' default. The special value -1 picks a number of threads based on the number
#' of cores in the system).
#'
#' @return A copy of the input data frame `df` with tje simulated values added
#' as a new column.
#'
#' @examples
#' # Assuming `fit` is a fitted GeoWarp model:
#' simulated_df <- geowarp_simulate(fit)
#'
#' @export
geowarp_simulate <- function(
  fit,
  df = fit$observed_df,
  model = fit$model,
  parameters = fit$parameters,
  n_samples = 1,
  nugget = TRUE,
  vecchia = 'auto',
  vecchia_n_parents = 50,
  parent_structure = vecchia_parent_structure(
    df,
    model,
    vecchia_n_parents
  ),
  threads = getOption('geowarp.threads')
) {
  .local_tbb_threads(threads)

  if (vecchia == 'auto') {
    vecchia <- nrow(df) > 1000
  }
  stan_data <- geowarp_stan_data(df, model)

  df[[model$variable]] <- if (!vecchia) {
    .geowarp_simulate_exact(
      stan_data,
      model,
      parameters,
      n_samples,
      nugget
    )
  } else {
    .geowarp_simulate_vecchia(
      stan_data,
      model,
      parameters,
      parent_structure,
      n_samples,
      nugget,
      threads
    )
  }
  df
}

.geowarp_simulate_exact <- function(
  stan_data,
  model,
  parameters,
  n_samples,
  nugget
) {
  Sigma <- .covariance_matrix_internal(
    stan_data$x,
    stan_data$X_deviation_fixed,
    stan_data$X_deviation_random,
    nugget = rep(nugget, nrow(stan_data$x)),
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

.geowarp_simulate_vecchia <- function(
  stan_data,
  model,
  parameters,
  parent_structure,
  n_samples,
  nugget,
  threads
) {
  U <- .create_vecchia_U(
    stan_data$x[parent_structure$observed_ordering, , drop = FALSE],
    stan_data$X_deviation_fixed[parent_structure$observed_ordering, , drop = FALSE],
    stan_data$X_deviation_random[parent_structure$observed_ordering, , drop = FALSE],
    parents = parent_structure$observed_parents,
    model = model,
    parameters = parameters,
    nugget = rep(nugget, nrow(stan_data$x)),
    threads = threads
  )

  output <- as.vector(cbind(
    stan_data$X_mean_fixed,
    stan_data$X_mean_random
  )[parent_structure$observed_ordering, , drop = FALSE] %*% parameters$alpha_beta) + solve(
    t(U),
    matrix(rnorm(n_samples * nrow(U)), ncol = n_samples)
  )
  as.matrix(output[Matrix::invPerm(parent_structure$observed_ordering), ])
}
