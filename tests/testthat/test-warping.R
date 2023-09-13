## bernstein_warping_design_matrix

test_that('bernstein_warping_design_matrix creates correct design matrix', {
  # Generate a simple example
  x <- seq(0, 5, length.out = 10)
  order <- 4
  domain <- c(0, 5)
  result <- bernstein_warping_design_matrix(x, order, domain)

  # Check dimensions
  expect_equal(nrow(result), length(x))
  expect_equal(ncol(result), order)

  # Check if the values are within expected range [0, 1]
  expect_true(all(result >= 0 & result <= 1))

  # Check first row of the matrix
  expect_equal(result[1, ], c(0, 0, 0, 0))
  # Check last row of the matrix
  expect_equal(result[10, ], c(0, 0, 0, 1))
  # Check middle row of the matrix
  x_prime <- x / domain[2]
  expect_equal(
    result[5, ],
    c(
      choose(order, 1) * x_prime[5] ^ 1 * (1 - x_prime[5]) ^ (order - 1),
      choose(order, 2) * x_prime[5] ^ 2 * (1 - x_prime[5]) ^ (order - 2),
      choose(order, 3) * x_prime[5] ^ 3 * (1 - x_prime[5]) ^ (order - 3),
      choose(order, 4) * x_prime[5] ^ 4 * (1 - x_prime[5]) ^ (order - 4)
    )
  )
})

test_that('bernstein_warping_design_matrix stops for invalid x values', {
  x <- c(-1, 0, 1)  # -1 is outside the domain [0, 1]
  order <- 2
  domain <- c(0, 1)

  expect_error(
    bernstein_warping_design_matrix(x, order, domain),
    'all\\(x_prime >= 0 & x_prime <= 1\\) is not TRUE'
  )
})

## warped_coordinates

make_warped_test_model <- function(deviation_model) {
  model <- geowarp_model(
    'value',
    c('x', 'y'),
    'z',
    list(c(0, 10), c(0, 10)),
    c(0, 10),
    deviation_model = deviation_model
  )
}
warped_test_input_data <- data.frame(
  x = c(0, 10),
  y = c(0, 10),
  z = c(0, 10)
)
warped_test_models <- list(
  linear_warpings = list(
    model = make_warped_test_model(
      geowarp_deviation_model(
        axial_warping_units = rep(list(geowarp_linear_awu()), 3),
        geometric_warping_unit = NULL
      )
    ),
    parameters = list(
      gamma_deviation_horizontal = c(1, 2),
      gamma_deviation_vertical = c(3)
    ),
    expected_second_row = c(10, 20, 30)
  ),
  full = list(
    model = make_warped_test_model(geowarp_deviation_model(
      axial_warping_units = list(
        geowarp_linear_awu(),
        geowarp_linear_awu(),
        geowarp_bernstein_awu()
      )
    )),
    parameters = list(
      gamma_deviation_horizontal = c(1, 2),
      gamma_deviation_vertical = c(rep(1, 9), 0.5, rep(1, 10)),
      # Technically not a valid L matrix
      L_deviation = rbind(
        c(1, 0, 0),
        c(0.75, 1, 0),
        c(0, 0, 1)
      )
    ),
    expected_second_row = c(25, 20, 19.5)
  ),
  white = list(
    model = make_warped_test_model(geowarp_white_deviation_model()),
    parameters = NULL,
    expected_second_row = c(10, 10, 10)
  ),
  vertical_only = list(
    model = make_warped_test_model(
      geowarp_vertical_only_deviation_model()
    ),
    parameters = list(
      gamma_deviation_vertical = c(rep(1, 9), 0.5, rep(1, 10))
    ),
    expected_second_row = c(10, 10, 19.5)
  )
)

test_that('warped_coordinates works with all cases', {
  for (test_case in warped_test_models) {
    result <- warped_coordinates(
      df = warped_test_input_data,
      model = test_case$model,
      parameters = test_case$parameters
    )
    expect_equal(result[1, ], c(0, 0, 0))
    expect_equal(result[2, ], test_case$expected_second_row)
  }
})

## unwarped_coordinates

test_that('unwarped_coordinates works with all cases', {
  for (test_case in warped_test_models) {
    warped_x <- warped_coordinates(
      df = warped_test_input_data,
      model = test_case$model,
      parameters = test_case$parameters
    )
    result <- unwarped_coordinates(
      warped_x = warped_x,
      model = test_case$model,
      parameters = test_case$parameters
    )
    expect_equal(result[1, ], c(0, 0, 0))
    expect_equal(result[2, ], c(10, 10, 10))
  }
})
