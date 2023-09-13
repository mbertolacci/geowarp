parent_knn_basic_checks <- function(x, result, n_parents) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = length(x))
  }

  expect_equal(ncol(result), n_parents + 1)
  expect_equal(nrow(result), nrow(x))
  expect_equal(result[, 1], seq_len(nrow(x)))

  # Correctness check
  expect_equal(result[1, ], c(1, rep(NA, n_parents)))
  expect_equal(result[2, ], c(2, 1, rep(NA, n_parents - 1)))
}

test_that('parent_knn returns correct nearest neighbours', {
  x <- matrix(rnorm(50), ncol = 2)
  n_parents <- 5
  result <- parent_knn(x, n_parents)
  parent_knn_basic_checks(x, result, n_parents)
})

test_that('parent_knn works with a vector', {
  x <- rnorm(100)
  n_parents <- 5
  result <- parent_knn(x, n_parents)
  parent_knn_basic_checks(x, result, n_parents)
})

test_that('parent_knn works when there are zero or one locations', {
  x <- matrix(rnorm(2), ncol = 2)
  n_parents <- 5
  result <- parent_knn(x, n_parents)
  expect_equal(result, rbind(c(1, rep(NA, n_parents))))

  x <- matrix(NA, ncol = 2, nrow = 0)
  n_parents <- 5
  result <- parent_knn(x, n_parents)
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), n_parents + 1)
})

test_that('parent_knn_across_groups returns correct nearest neighbours from different groups', {
  x <- matrix(rnorm(100), nrow = 50)
  n_parents <- 5
  groups <- c(1, 1, sample(1 : 5, nrow(x) - 2, replace = TRUE))
  result <- parent_knn_across_groups(x, groups, n_parents)

  expect_equal(ncol(result), n_parents + 1)
  expect_equal(nrow(result), nrow(x))
  expect_equal(result[, 1], seq_len(nrow(x)))

  expect_equal(result[1, ], c(1, rep(NA, n_parents)))
  # Second element is in the same group as the first, so first is not its parent
  expect_equal(result[2, ], c(2, rep(NA, n_parents)))
})

test_that('parent_knn_across_groups works with a vector', {
  x <- rnorm(100)
  n_parents <- 5
  groups <- c(1, 1, sample(1 : 5, length(x) - 2, replace = TRUE))
  result <- parent_knn_across_groups(x, groups, n_parents)

  expect_equal(ncol(result), n_parents + 1)
  expect_equal(nrow(result), length(x))
  expect_equal(result[, 1], seq_along(x))

  expect_equal(result[1, ], c(1, rep(NA, n_parents)))
  # Second element is in the same group as the first, so first is not its parent
  expect_equal(result[2, ], c(2, rep(NA, n_parents)))
})

test_that('parent_knn_across_groups fails when the number of groups doesn\'t match', {
  x <- matrix(rnorm(100), nrow = 50)
  n_parents <- 5
  groups <- rep(1, 10)
  expect_error(parent_knn_across_groups(x, groups, n_parents))
})
