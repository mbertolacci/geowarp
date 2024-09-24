#' Calculate estimated depth-dependent mean and variance
#'
#' In the GeoWarp model, the process mean and variance vary with depth. These
#' functions calculate the estimated mean and variance at the given locations.
#'
#' @param fit A GeoWarp fit object form \code{\link{geowarp_optimise}}. Can be
#' omitted if both `model` and `parameters` are given.
#' @param model GeoWarp model object (default uses the model from the `fit`
#' object).
#' @param parameters GeoWarp parameters (default uses the parameters from the
#' `fit` object).
#' @param df A data frame containing the observed data (default uses the
#' observation locations used in \code{fit}).
#' @param nugget Logical indicating whether to include the nugget variance in
#' the output (default is TRUE).
#' @return A vector containing the estimated mean or marginal variance at each
#' location in `df`.
#'
#' @examples
#' # Assuming `fit` is a fitted GeoWarp model:
#' mean_values <- mean_profile(fit)
#' @export
mean_profile <- function(
  fit,
  df = fit$observed_df,
  model = fit$model,
  parameters = fit$parameters
) {
  stan_data <- geowarp_stan_data(df, model)

  as.vector(cbind(
    stan_data$X_mean_fixed,
    stan_data$X_mean_random
  ) %*% parameters$alpha_beta)
}

#' @describeIn mean_profile Output the marginal variance at each depth
#' @export
marginal_variance_profile <- function(
  fit,
  df = fit$observed_df,
  model = fit$model,
  parameters = fit$parameters,
  nugget = TRUE
) {
  stan_data <- geowarp_stan_data(df, model)

  as.vector(exp(
    stan_data$X_deviation_fixed %*% parameters$eta_deviation
    + stan_data$X_deviation_random %*% parameters$zeta_deviation
  ) + nugget * parameters$sigma_squared_nugget)
}

#' Predict from a Fitted GeoWarp Model
#'
#' This function generates predictions using a fitted GeoWarp model. The
#' prediction can use the Vecchia approximation when the number of observed and
#' predicted data are large. The prediction distribution is multivariate
#' Gaussian, so it is specified by its mean and covariance
#'
#' @param fit GeoWarp model fit produced using \code{\link{geowarp_optimise}}.
#' @param prediction_df Data frame for locations where predictions are to be
#' made. Defaults to the observed data frame in the model fit.
#' @param observed_df Data frame containing the observed data. Defaults to the
# observed data frame in the model fit.
#' @param model GeoWarp model object (default uses the model from the fit).
#' @param parameters GeoWarp parameters (default uses the parameters from the
#' fit).
#' @param include_mean Logical, indicating whether to include the mean in the
#' output.
#' @param include_precision_V Logical, indicating whether to include the V
#' matrix in the output (explained below).
#' @param include_samples Logical, indicating whether to include the predictive
#' samples in the output.
#' @param nugget Logical, indicating whether to include a nugget effect in the
#' prediction.
#' @param n_samples Number of predictive samples to generate.
#' @param vecchia Whether to use Vecchia's approximation (`'auto'`, `TRUE`,
#' `FALSE`). If `'auto'`, the Vecchia approximation is used when the number of
#' observed and prediction locations together exceed 1000.
#' @param n_parents Number of parents to consider for Vecchia's approximation.
#' @param parent_structure Parent structure for Vecchia's approximation,
#' produced using \code{\link{vecchia_parent_structure}}.
#' @param threads Number of threads to use for parallelisation, taken from
#' `getOption('geowarp.threads')` by default, which is itself set to 1 by
#' default. The special value -1 picks a number of threads based on the number
#' of cores in the system). Only applies if using the Vecchia approximation.
#' @param ... Ignored.
#'
#' @return A list containing the elements of the predictive distribution
#' (depending on the options):
#' \itemize{
#' \item \code{mean}: The mean of the predictive distribution.
#' \item \code{precision_V}: A (possibly sparse) matrix V such that VV' is the
#' reordered precision matrix of the predictive distribution.
#' \item \code{ordering}: A vector of integers indicating the ordering of the
#' entries in the matrix V; see Details below.
#' \item \code{samples}: A matrix of predictive samples.
#' }
#'
#' @section Details:
#' If `prediction` contains the output of this function, the precision matrix
#' of the prediction can be recovered by calculating
#' `tcrossprod(prediction$precision_V)[prediction$ordering, prediction$ordering]`.
#' However, this may be of limited use because, while the precision matrix may
#' be sparse, the covariance may not be. To calculate variances and other
#' quantities it is better to draw samples from the prediction distribution.
#' This function can do that, or you can do it manually as shown in the example
#' section.
#'
#' @examples
#' # Assuming `fit` is a fitted GeoWarp model and `new_data` is a new data frame
#' # for prediction:
#' prediction_dist <- predict(
#'   fit,
#'   prediction_df = new_data,
#'   include_precision_V = TRUE,
#'   include_samples = TRUE
#' )
#' # Samples are now in `prediction_dist$samples`, but you can also generate
#' # them manually
#' y_tilde <- as.matrix(solve(
#'   t(prediction_dist$precision_V),
#'   z
#' ))
#' manual_samples <- prediction_dist$mean + y_tilde[prediction_dist$ordering, ]
#'
#' @seealso
#' \code{\link{geowarp_optimise}}
#' \code{\link{vecchia_parent_structure}}
#' @export
geowarp_predict <- function(
  fit,
  prediction_df = fit$observed_df,
  observed_df = fit$observed_df,
  parameters = fit$parameters,
  model = fit$model,
  include_mean = TRUE,
  include_precision_V = FALSE,
  include_samples = FALSE,
  nugget = TRUE,
  n_samples = 500,
  vecchia = 'auto',
  n_parents = 50,
  parent_structure = vecchia_parent_structure(
    observed_df,
    model,
    n_parents,
    prediction_df = prediction_df
  ),
  threads = getOption('geowarp.threads'),
  ...
) {
  .local_tbb_threads(threads)

  if (vecchia == 'auto') {
    vecchia <- (nrow(observed_df) + nrow(prediction_df)) > 1000
  }

  if (model$deviation_model$name == 'white') {
    prediction_distribution <- .prediction_distribution_white(
      observed_df,
      prediction_df,
      model,
      parameters
    )
  } else if (!vecchia) {
    prediction_distribution <- .prediction_distribution_exact(
      observed_df,
      prediction_df,
      model,
      parameters,
      nugget
    )
  } else {
    prediction_distribution <- .prediction_distribution_vecchia(
      observed_df,
      prediction_df,
      model,
      parameters,
      parent_structure,
      nugget,
      threads
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
    z <- matrix(rnorm(nrow(prediction_df) * n_samples), nrow = nrow(prediction_df))
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

#' @describeIn geowarp_predict Predict using a fit object
#' @export
predict.geowarp_fit <- function(object, ...) {
  geowarp_predict(object, ...)
}

.prediction_distribution_exact <- function(
  observed_df,
  prediction_df,
  model,
  parameters,
  nugget
) {
  observed_stan_data <- geowarp_stan_data(observed_df, model)
  prediction_stan_data <- geowarp_stan_data(prediction_df, model)

  rbind_matrices <- function(name) {
    rbind(observed_stan_data[[name]], prediction_stan_data[[name]])
  }
  Sigma_all <- .covariance_matrix_internal(
    rbind_matrices('x'),
    rbind_matrices('X_deviation_fixed'),
    rbind_matrices('X_deviation_random'),
    model = model,
    parameters = parameters,
    nugget = c(
      rep(TRUE, nrow(observed_stan_data$x)),
      rep(nugget, nrow(prediction_stan_data$x))
    )
  )

  n_observed <- nrow(observed_stan_data$x)
  n_prediction <- nrow(prediction_stan_data$x)
  observed_indices <- seq_len(n_observed)
  prediction_indices <- n_observed + seq_len(n_prediction)

  Sigma_observed <- Sigma_all[observed_indices, observed_indices]
  Sigma_prediction <- Sigma_all[prediction_indices, prediction_indices]
  Sigma_cross <- Sigma_all[observed_indices, prediction_indices]

  mu_observed <- as.vector(cbind(
    observed_stan_data$X_mean_fixed,
    observed_stan_data$X_mean_random
  ) %*% parameters$alpha_beta)
  mu_prediction <- as.vector(cbind(
    prediction_stan_data$X_mean_fixed,
    prediction_stan_data$X_mean_random
  ) %*% parameters$alpha_beta)

  chol_Sigma_observed <- chol(Sigma_observed)
  mu_prediction <- mu_prediction + as.vector(crossprod(Sigma_cross, .chol_solve(
    chol_Sigma_observed,
    observed_stan_data$y - mu_observed
  )))
  Sigma_prediction <- (
    Sigma_prediction
    - crossprod(Sigma_cross, .chol_solve(chol_Sigma_observed, Sigma_cross))
  )
  Q_prediction <- chol2inv(chol(Sigma_prediction))
  rev_matrix <- function(x) {
    x[rev(seq_len(nrow(x))), rev(seq_len(ncol(x)))]
  }
  list(
    mean = mu_prediction,
    precision_V = rev_matrix(t(chol(rev_matrix(Q_prediction)))),
    ordering = seq_len(n_prediction)
  )
}

.prediction_distribution_vecchia <- function(
  observed_df,
  prediction_df,
  model,
  parameters,
  parent_structure,
  nugget,
  threads
) {
  observed_stan_data <- geowarp_stan_data(observed_df, model)
  prediction_stan_data <- geowarp_stan_data(prediction_df, model)

  rbind_matrices <- function(name) {
    rbind(
      observed_stan_data[[name]][parent_structure$observed_ordering, , drop = FALSE],
      prediction_stan_data[[name]][parent_structure$prediction_ordering, , drop = FALSE]
    )
  }
  x_all <- rbind_matrices('x')
  X_mean_all <- cbind(
    rbind_matrices('X_mean_fixed'),
    rbind_matrices('X_mean_random')
  )
  X_deviation_fixed_all <- rbind_matrices('X_deviation_fixed')
  X_deviation_random_all <- rbind_matrices('X_deviation_random')

  n_observed <- nrow(observed_stan_data$x)
  n_prediction <- nrow(prediction_stan_data$x)
  observed_indices <- seq_len(n_observed)
  prediction_indices <- n_observed + seq_len(n_prediction)

  parents <- rbind(
    parent_structure$observed_parents,
    .deduplicate_parents(cbind(
      nrow(parent_structure$observed_parents) + parent_structure$prediction_within_parents,
      parent_structure$prediction_between_parents
    ))
  )
  U <- .create_vecchia_U(
    x_all,
    X_deviation_fixed_all,
    X_deviation_random_all,
    parents = parents,
    model = model,
    parameters = parameters,
    nugget = c(
      rep(TRUE, nrow(observed_stan_data$x)),
      rep(nugget, nrow(prediction_stan_data$x))
    ),
    threads = threads
  )
  V <- U[prediction_indices, prediction_indices]

  mu_hat <- as.vector(X_mean_all %*% parameters$alpha_beta)
  prediction_mean <- (
    mu_hat[prediction_indices]
    - as.vector(solve(t(V), solve(V,
      U[prediction_indices, ] %*%
      crossprod(
        U[observed_indices, ],
        observed_stan_data$y[parent_structure$observed_ordering]
        - mu_hat[observed_indices]
      )
    )))
  )
  reverse_ordering <- invPerm(parent_structure$prediction_ordering)
  list(
    mean = prediction_mean[reverse_ordering],
    precision_V = V,
    ordering = reverse_ordering
  )
}

.prediction_distribution_white <- function(
  observed_df,
  prediction_df,
  model,
  parameters
) {
  prediction_stan_data <- geowarp_stan_data(prediction_df, model)

  prediction_mean <- as.vector(cbind(
    prediction_stan_data$X_mean_fixed,
    prediction_stan_data$X_mean_random
  ) %*% parameters$alpha_beta)
  prediction_sd <- exp(as.vector(0.5 * (
    prediction_stan_data$X_deviation_fixed %*% parameters$eta_deviation
    + prediction_stan_data$X_deviation_random %*% parameters$zeta_deviation
  )))

  precision_V <- Diagonal(x = 1 / prediction_sd)
  list(
    mean = prediction_mean,
    precision_V = precision_V,
    ordering = seq_len(prediction_stan_data$N)
  )
}

#' Give the posterior distribution of alpha beta
#'
#' This function gives the posterior distribution of the alpha beta parameters
#' in a fitted GeoWarp model.
#'
#' @param fit A GeoWarp fit object from \code{\link{geowarp_optimise}}.
#' @param df A data frame containing the observed data (default uses the
#' observation locations used in \code{fit}).
#' @param parameters GeoWarp parameters (default uses the parameters from the
#' \code{fit} object).
#' @param model GeoWarp model object (default uses the model from the \code{fit}
#' object).
#' @param include_mean Logical, indicating whether to include the mean in the
#' output.
#' @param include_precision Logical, indicating whether to include the precision
#' matrix in the output (explained below).
#' @param include_samples Logical, indicating whether to include the predictive
#' samples in the output.
#' @param n_samples Number of predictive samples to generate.
#' @param n_parents Number of parents to consider for Vecchia's approximation.
#' @param parent_structure Parent structure for Vecchia's approximation,
#' produced using \code{\link{vecchia_parent_structure}}.
#' @param threads Number of threads to use for parallelisation, taken from
#' \code{getOption('geowarp.threads')} by default, which is itself set to 1 by
#' default. The special value -1 picks a number of threads based on the number
#' of cores in the system).
#'
#' @return A vector containing the sampled alpha beta parameters.
#' @export
predict_alpha_beta <- function(
  fit,
  df = fit$observed_df,
  parameters = fit$parameters,
  model = fit$model,
  include_mean = TRUE,
  include_precision = FALSE,
  include_samples = FALSE,
  n_samples = 1L,
  n_parents = 50,
  parent_structure = vecchia_parent_structure(
    df,
    model,
    n_parents
  ),
  threads = getOption('geowarp.threads')
) {
  stan_data <- geowarp_stan_data(df, model)
  stan_data$y <- stan_data$y[parent_structure$observed_ordering]
  for (name in c(
    'x', 'X_mean_fixed', 'X_mean_random', 'X_deviation_warping',
    'X_deviation_fixed', 'X_deviation_random'
  )) {
    stan_data[[name]] <- stan_data[[name]][
      parent_structure$observed_ordering,
      ,
      drop = FALSE
    ]
  }

  X_mean_all <- cbind(stan_data$X_mean_fixed, stan_data$X_mean_random)

  U <- .create_vecchia_U(
    x = stan_data$x,
    X_deviation_fixed = stan_data$X_deviation_fixed,
    X_deviation_random = stan_data$X_deviation_random,
    parents = parent_structure$observed_parents,
    nugget = rep(TRUE, nrow(stan_data$x)),
    model = model,
    parameters = parameters,
    threads = threads
  )
  alpha_beta_precision <- (
    as.matrix(crossprod(
      X_mean_all,
      as.matrix(
        U %*% as.matrix(crossprod(U, X_mean_all))
      )
    ))
    + as.matrix(Matrix::bdiag(
      stan_data$alpha_precision,
      .rw1d_precision(
        ncol(stan_data$X_mean_random),
        parameters$tau_squared_mean_random
      )
    ))
  )
  if (include_mean || include_samples) {
    alpha_beta_hat <- as.vector(solve(
      alpha_beta_precision,
      crossprod(X_mean_all, as.matrix(U %*% crossprod(U, stan_data$y)))
      + c(
        stan_data$alpha_precision %*% stan_data$alpha_mean,
        rep(0, ncol(stan_data$X_mean_random))
      )
    ))
  }

  output <- NULL
  if (include_mean) {
    output$mean <- alpha_beta_hat
  }
  if (include_precision) {
    output$precision <- alpha_beta_precision
  }
  if (include_samples) {
    output$samples <- t(alpha_beta_hat + backsolve(
      chol(alpha_beta_precision),
      matrix(
        rnorm(n_samples * nrow(alpha_beta_precision)),
        nrow = nrow(alpha_beta_precision),
        ncol = n_samples
      )
    ))
    if (n_samples == 1) {
      output$samples <- output$samples[1, ]
    }
  }
  output
}
