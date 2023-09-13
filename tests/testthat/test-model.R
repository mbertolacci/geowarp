test_that('default geowarp_model call works', {
  model <- geowarp_model(
    variable = 'value',
    horizontal_coordinates = 'x',
    vertical_coordinate = 'z',
    horizontal_domains = c(0, 10),
    vertical_domain = c(0, 10),
  )
  expect_equal(model$horizontal_domains, list(c(0, 10)))
})

test_that('inputs are validated correctly', {
  expect_error({
    geowarp_model(
      horizontal_coordinates = c('x', 'y'),
      horizontal_domains = list(c(0, 10), c(10, 0))
    )
  })
  expect_error({
    geowarp_model(
      horizontal_coordinates = c('x', 'y'),
      horizontal_domains = list(c(0, 10))
    )
  })
  expect_error({
    geowarp_model(
      horizontal_coordinates = c('x', 'y'),
      horizontal_domains = list(c(0, 10), c(0, 10)),
      vertical_domain = c(10, 0)
    )
  })
})

test_that('geowarp_mean_model can take a formula or a string', {
  model1 <- geowarp_mean_model(~ depth)
  expect_true(inherits(model1$fixed_formula, 'formula'))
  model2 <- geowarp_mean_model('depth')
  expect_true(inherits(model2$fixed_formula, 'formula'))
})

test_that('geowarp_deviation_model checks warpings correctly', {
  # Valid; no error
  model <- geowarp_deviation_model(axial_warping_units = list(
    geowarp_linear_awu(),
    geowarp_bernstein_awu()
  ))
  expect_equal(model$name, 'full')
  expect_error({
    # Horizontal warpings must be linear
    geowarp_deviation_model(axial_warping_units = list(
      geowarp_bernstein_awu(),
      geowarp_bernstein_awu()
    ))
  })
})

test_that('geowarp_white_deviation_model and geowarp_vertical_only_deviation_model run cleanly', {
  expect_equal(
    geowarp_white_deviation_model()$name,
    'white'
  )
  expect_equal(
    geowarp_vertical_only_deviation_model()$name,
    'vertical_only'
  )
})

test_that('geowarp_variance_model runs cleanly', {
  expect_equal(
    geowarp_variance_model()$fixed_effect_mean,
    0
  )
})

test_that('warping units run cleanly', {
  expect_equal(geowarp_linear_awu()$name, 'linear_awu')
  expect_equal(geowarp_bernstein_awu()$name, 'bernstein_awu')
  expect_equal(geowarp_geometric_warping_unit()$name, 'geometric_warping_unit')
})

test_that('model_coordinates extracts the right columns', {
  model <- geowarp_model(
    variable = 'value',
    horizontal_coordinates = c('x', 'y'),
    horizontal_domains = list(c(0, 10), c(0, 10)),
    vertical_coordinate = 'z',
    vertical_domain = c(0, 10)
  )
  result <- model_coordinates(data.frame(
    x = rep(1, 2),
    y = rep(2, 2),
    z = rep(3, 2)
  ), model)
  expect_equal(result, cbind(rep(1, 2), rep(2, 2), rep(3, 2)))
})
