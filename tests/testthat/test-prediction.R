check_mean_profile_and_marginal_variance_profile <- function(fn) {
  model <- geowarp_model(
    variable = 'value',
    horizontal_coordinates = 'x',
    vertical_coordinate = 'z',
    horizontal_domains = c(0, 10),
    vertical_domain = c(0, 10)
  )
  target_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  result <- fn(
    df = target_df,
    model = model,
    parameters = list(
      alpha_beta = rep(0, 105),
      eta_deviation = 0,
      zeta_deviation = rep(0, 13),
      gamma_deviation_horizontal = 1,
      gamma_deviation_vertical = rep(0.1, 20),
      L_deviation = diag(2),
      sigma_squared_nugget = 1
    )
  )
  expect_equal(length(result), nrow(target_df))
  expect_type(result, 'double')
}

## mean_profile

test_that('mean_profile works', {
  check_mean_profile_and_marginal_variance_profile(mean_profile)
})

## marginal_variance_profile

test_that('marginal_variance_profile works', {
  check_mean_profile_and_marginal_variance_profile(marginal_variance_profile)
})

## geowarp_predict

full_model <- geowarp_model(
  variable = 'value',
  horizontal_coordinates = 'x',
  vertical_coordinate = 'z',
  horizontal_domains = c(0, 10),
  vertical_domain = c(0, 10)
)
observed_df <- data.frame(
  x = rep(1 : 5, each = 5),
  z = rep(1 : 5, 5),
  value = rnorm(25)
)
parameters <- list(
  alpha_beta = rep(0, 105),
  eta_deviation = 0,
  zeta_deviation = rep(0, 13),
  gamma_deviation_horizontal = 1,
  gamma_deviation_vertical = rep(0.1, 20),
  L_deviation = diag(2),
  sigma_squared_nugget = 1
)

run_predict <- function(model = full_model, ...) {
  geowarp_predict(
    prediction_df = observed_df,
    observed_df = observed_df,
    parameters = parameters,
    model = model,
    ...
  )
}

test_that('prediction output flags work', {
  result <- run_predict(include_mean = TRUE, include_precision_V = FALSE, include_samples = FALSE)
  expect_equal(length(result$mean), nrow(observed_df))
  expect_null(result$precision_V)
  expect_null(result$samples)

  result <- run_predict(include_mean = FALSE, include_precision_V = TRUE, include_samples = FALSE)
  expect_null(result$mean)
  expect_equal(dim(result$precision_V), rep(nrow(observed_df), 2))

  result <- run_predict(include_mean = FALSE, include_precision_V = FALSE, include_samples = TRUE, n_samples = 10)
  expect_equal(dim(result$samples), c(nrow(observed_df), 10))
})

test_that('prediction with a white model works', {
  result <- run_predict(
    model = geowarp_model(
      variable = 'value',
      horizontal_coordinates = 'x',
      vertical_coordinate = 'z',
      horizontal_domains = c(0, 10),
      vertical_domain = c(0, 10),
      deviation_model = geowarp_white_deviation_model()
    ),
    include_mean = TRUE,
    include_precision_V = TRUE,
    include_samples = TRUE,
    n_samples = 10
  )
  expect_equal(length(result$mean), nrow(observed_df))
  expect_equal(dim(result$precision_V), rep(nrow(observed_df), 2))
  expect_equal(dim(result$samples), c(nrow(observed_df), 10))
})

test_that('prediction with a vertical only model works', {
  result <- run_predict(
    model = geowarp_model(
      variable = 'value',
      horizontal_coordinates = 'x',
      vertical_coordinate = 'z',
      horizontal_domains = c(0, 10),
      vertical_domain = c(0, 10),
      deviation_model = geowarp_vertical_only_deviation_model()
    ),
    include_mean = TRUE,
    include_precision_V = TRUE,
    include_samples = TRUE,
    n_samples = 10
  )
  expect_equal(length(result$mean), nrow(observed_df))
  expect_equal(dim(result$precision_V), rep(nrow(observed_df), 2))
  expect_equal(dim(result$samples), c(nrow(observed_df), 10))
})

test_that('prediction with the vecchia approximation works', {
  result <- run_predict(
    vecchia = TRUE,
    include_mean = TRUE,
    include_precision_V = TRUE,
    include_samples = TRUE,
    n_samples = 10
  )
  expect_equal(length(result$mean), nrow(observed_df))
  expect_equal(dim(result$precision_V), rep(nrow(observed_df), 2))
  expect_s4_class(result$precision_V, 'sparseMatrix')
  expect_equal(dim(result$samples), c(nrow(observed_df), 10))
})

test_that('prediction with the vecchia approximation works with threads > 1', {
  result <- run_predict(
    vecchia = TRUE,
    include_mean = TRUE,
    include_precision_V = TRUE,
    include_samples = TRUE,
    n_samples = 10,
    threads = 2
  )
  expect_equal(length(result$mean), nrow(observed_df))
  expect_equal(dim(result$precision_V), rep(nrow(observed_df), 2))
  expect_s4_class(result$precision_V, 'sparseMatrix')
  expect_equal(dim(result$samples), c(nrow(observed_df), 10))
})

test_that('prediction with the vecchia approximation and a vertical only model works', {
  result <- run_predict(
    model = geowarp_model(
      variable = 'value',
      horizontal_coordinates = 'x',
      vertical_coordinate = 'z',
      horizontal_domains = c(0, 10),
      vertical_domain = c(0, 10),
      deviation_model = geowarp_vertical_only_deviation_model()
    ),
    vecchia = TRUE,
    include_mean = TRUE,
    include_precision_V = TRUE,
    include_samples = TRUE,
    n_samples = 10
  )
  expect_equal(length(result$mean), nrow(observed_df))
  expect_equal(dim(result$precision_V), rep(nrow(observed_df), 2))
  expect_s4_class(result$precision_V, 'sparseMatrix')
  expect_equal(dim(result$samples), c(nrow(observed_df), 10))
})

## predict.geowarp_fit

test_that('predict.geowarp_fit works', {
  fit_shim <- list(
    parameters = parameters,
    model = full_model,
    observed_df = observed_df
  )
  class(fit_shim) <- 'geowarp_fit'
  result <- predict(fit_shim)
  expect_equal(length(result$mean), nrow(observed_df))
})
