test_that('gamma_quantile_prior returns correct shape and rate parameters', {
  result <- gamma_quantile_prior(q_lower = 2, q_upper = 10)
  expect_type(result, 'list')
  expect_true(all(names(result) %in% c('shape', 'rate')))

  calculated_quantiles <- qgamma(c(0.05, 0.95), shape = result$shape, rate = result$rate)
  expect_equal(calculated_quantiles[1], 2, tolerance = 1e-5)
  expect_equal(calculated_quantiles[2], 10, tolerance = 1e-5)
})

test_that('inverse_gamma_quantile_prior returns correct shape and rate parameters', {
  result <- inverse_gamma_quantile_prior(q_lower = 2, q_upper = 10)
  expect_type(result, 'list')
  expect_true(all(names(result) %in% c('shape', 'rate')))

  calculated_quantiles <- qgamma(c(0.95, 0.05), shape = result$shape, rate = result$rate)
  expect_equal(1 / calculated_quantiles[1], 2, tolerance = 1e-5)
  expect_equal(1 / calculated_quantiles[2], 10, tolerance = 1e-5)
})
