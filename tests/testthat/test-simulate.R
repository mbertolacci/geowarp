## geowarp_simulate

test_that('simulation works', {
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

  call_simulate <- function(...) {
    geowarp_simulate(
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
      ),
      ...
    )
  }

  base_checks <- function(result) {
    expect_true(is.data.frame(result))
    expect_equal(colnames(result), c('x', 'z', 'value'))
    expect_equal(result$x, target_df$x)
    expect_equal(result$z, target_df$z)
  }

  result <- call_simulate(n_sample = 1, vecchia = 'auto')
  base_checks(result)
  expect_type(result$value, 'double')

  result <- call_simulate(n_sample = 1, vecchia = FALSE)
  base_checks(result)
  expect_type(result$value, 'double')

  result <- call_simulate(n_sample = 5, vecchia = FALSE)
  base_checks(result)
  expect_equal(dim(result$value), c(nrow(target_df), 5))

  result <- call_simulate(n_sample = 1, vecchia = TRUE)
  base_checks(result)
  expect_type(result$value, "double")
  
  result <- call_simulate(n_sample = 5, vecchia = TRUE)
  base_checks(result)
  expect_equal(dim(result$value), c(nrow(target_df), 5))

  result <- call_simulate(n_sample = 1, vecchia = TRUE, threads = 2L)
  base_checks(result)
})
