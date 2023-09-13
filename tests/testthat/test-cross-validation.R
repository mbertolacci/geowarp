simple_model <- geowarp_model(
  variable = 'value',
  horizontal_coordinates = 'x',
  vertical_coordinate = 'z',
  horizontal_domains = c(0, 10),
  vertical_domain = c(0, 10),
  mean_model = geowarp_mean_model(
    fixed_formula = ~ 1,
    vertical_basis_functions = FALSE
  ),
  deviation_model = geowarp_deviation_model(
    axial_warping_units = rep(list(geowarp_linear_awu()), 2),
    variance_model = geowarp_variance_model(
      vertical_basis_functions = FALSE
    )
  )
)

observed_df <- data.frame(
  group = sample(1 : 5, 30, replace = TRUE),
  z = runif(30, 0, 10),
  value = rnorm(30)
) %>%
  mutate(x = group + 0.1)

test_that('geowarp_cross_validation_fit with a simple model works', {
  fits <- geowarp_cross_validation_fit(
    observed_df,
    'group',
    simple_model
  )
  expect_equal(length(fits), 5)
})

test_that('geowarp_cross_validation_fit works with show_progress = TRUE', {
  fits <- geowarp_cross_validation_fit(
    observed_df,
    'group',
    simple_model,
    show_progress = TRUE
  )
  expect_equal(length(fits), 5)
})

test_that('geowarp_cross_validation_fit output works', {
  path <- tempdir()
  fits <- geowarp_cross_validation_fit(
    observed_df,
    'group',
    simple_model,
    group_output_pattern = file.path(path, 'fit_{group}.txt')
  )
  expect_true(all(file.exists(file.path(path, sprintf('fit_%d.txt', 1 : 5)))))
})
