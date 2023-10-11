simple_model <- geowarp_model(
  variable = 'value',
  horizontal_coordinates = 'x',
  vertical_coordinate = 'z',
  horizontal_domains = c(0, 10),
  vertical_domain = c(0, 10),
  mean_model = geowarp_mean_model(
    fixed_formula = ~ 1,
    vertical_basis_function_delta = 5
  ),
  deviation_model = geowarp_deviation_model(
    axial_warping_units = rep(list(geowarp_linear_awu()), 2),
    variance_model = geowarp_variance_model(
      vertical_basis_function_delta = 5
    )
  )
)

observed_df <- data.frame(
  x = runif(30, 0, 10),
  z = runif(30, 0, 10),
  value = rnorm(30)
)

test_that('geowarp_optimise works with simple model', {
  fit <- geowarp_optimise(observed_df, simple_model, best_of = 2)
  expect_equal(length(fit$stan_fit$attempts), 2)
  expect_equal(length(fit$parameters$tau_squared_mean_random), 1)
  expect_equal(length(fit$parameters$sigma_squared_nugget), 1)
  expect_equal(length(fit$parameters$eta_deviation), 1)
  expect_equal(length(fit$parameters$zeta_deviation), 5)
  expect_equal(length(fit$parameters$ell_deviation_random), 1)
  expect_equal(length(fit$parameters$tau_squared_deviation_random), 1)
  expect_equal(length(fit$parameters$gamma_deviation_horizontal), 1)
  expect_equal(length(fit$parameters$gamma_deviation_vertical), 1)
  expect_equal(dim(fit$parameters$L_deviation), c(2, 2))
  expect_equal(length(fit$parameters$alpha_beta), 6)
})

test_that('geowarp_optimise works with threading enabled', {
  fit <- geowarp_optimise(observed_df, simple_model, best_of = 2, threads = 2)
  expect_equal(length(fit$stan_fit$attempts), 2)
})

test_that('geowarp_optimise works with optimizing method', {
  fit <- geowarp_optimise(observed_df, simple_model, best_of = 2, method = 'optimizing')
  expect_equal(length(fit$stan_fit$attempts), 2)
})

test_that('geowarp_optimise works with fixed effect precision matrices', {
  simple_model$mean_model$fixed_effect_precision <- as.matrix(
    simple_model$mean_model$fixed_effect_precision
  )
  simple_model$deviation_model$variance_model$fixed_effect_precision <- as.matrix(
    simple_model$deviation_model$variance_model$fixed_effect_precision
  )
  fit <- geowarp_optimise(observed_df, simple_model, best_of = 2)
  expect_equal(length(fit$stan_fit$attempts), 2)
})

test_that('geowarp_optimise data validation works', {
  expect_error(
    geowarp_optimise(
      bind_rows(observed_df, observed_df),
      simple_model,
      best_of = 2
    ),
    'Duplicate coordinates are not supported'
  )
  expect_error(
    geowarp_optimise(
      observed_df %>% mutate(x = x + 100),
      simple_model,
      best_of = 2
    ),
    'Coordinate x is outside of the domain'
  )
})

test_that('geowarp_optimise works with vecchia approximation', {
  fit <- geowarp_optimise(
    observed_df,
    simple_model,
    n_parents = 5,
    vecchia = TRUE,
    best_of = 2
  )
  expect_equal(length(fit$stan_fit$attempts), 2)
  expect_equal(length(fit$parameters$tau_squared_mean_random), 1)
  expect_equal(length(fit$parameters$sigma_squared_nugget), 1)
  expect_equal(length(fit$parameters$eta_deviation), 1)
  expect_equal(length(fit$parameters$zeta_deviation), 5)
  expect_equal(length(fit$parameters$ell_deviation_random), 1)
  expect_equal(length(fit$parameters$tau_squared_deviation_random), 1)
  expect_equal(length(fit$parameters$gamma_deviation_horizontal), 1)
  expect_equal(length(fit$parameters$gamma_deviation_vertical), 1)
  expect_equal(dim(fit$parameters$L_deviation), c(2, 2))
  expect_equal(length(fit$parameters$alpha_beta), 6)
})

test_that('geowarp_optimise works with vertical only model', {
  vertical_only_model <- geowarp_model(
    variable = 'value',
    horizontal_coordinates = 'x',
    vertical_coordinate = 'z',
    horizontal_domains = c(0, 10),
    vertical_domain = c(0, 10),
    mean_model = geowarp_mean_model(
      fixed_formula = ~ 1,
      vertical_basis_functions = FALSE
    ),
    deviation_model = geowarp_vertical_only_deviation_model()
  )
  fit <- geowarp_optimise(observed_df, vertical_only_model, best_of = 2)
  expect_equal(length(fit$stan_fit$attempts), 2)
  expect_equal(length(fit$parameters$tau_squared_mean_random), 1)
  expect_equal(length(fit$parameters$sigma_squared_nugget), 1)
  expect_equal(length(fit$parameters$eta_deviation), 1)
  expect_equal(length(fit$parameters$zeta_deviation), 13)
  expect_equal(length(fit$parameters$ell_deviation_random), 1)
  expect_equal(length(fit$parameters$gamma_deviation_vertical), 20)
  expect_equal(length(fit$parameters$tau_squared_deviation_random), 1)
  expect_equal(length(fit$parameters$alpha_beta), 1)
})

test_that('geowarp_optimise works with white deviation model', {
  white_model <- geowarp_model(
    variable = 'value',
    horizontal_coordinates = 'x',
    vertical_coordinate = 'z',
    horizontal_domains = c(0, 10),
    vertical_domain = c(0, 10),
    mean_model = geowarp_mean_model(
      fixed_formula = ~ 1,
      vertical_basis_functions = FALSE
    ),
    deviation_model = geowarp_white_deviation_model()
  )
  fit <- geowarp_optimise(observed_df, white_model, best_of = 2)
  expect_equal(length(fit$stan_fit$attempts), 2)
  expect_equal(length(fit$parameters$tau_squared_mean_random), 1)
  expect_equal(length(fit$parameters$sigma_squared_nugget), 1)
  expect_equal(length(fit$parameters$eta_deviation), 1)
  expect_equal(length(fit$parameters$zeta_deviation), 13)
  expect_equal(length(fit$parameters$ell_deviation_random), 1)
  expect_equal(length(fit$parameters$tau_squared_deviation_random), 1)
  expect_equal(length(fit$parameters$alpha_beta), 1)
})
