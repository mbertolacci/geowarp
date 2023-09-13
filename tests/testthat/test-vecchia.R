make_model <- function(...) {
  geowarp_model(
    variable = 'value',
    horizontal_coordinates = 'x',
    vertical_coordinate = 'z',
    horizontal_domains = c(0, 10),
    vertical_domain = c(0, 10),
    ...
  )
}

full_model <- make_model()
vertical_only_model <- make_model(
  deviation_model = geowarp_vertical_only_deviation_model()
)

## vecchia_parent_structure
# Note: not attempting to check correctness here, just running the code through
# its paces

test_that('vecchia_parent_structure fails for white deviation model', {
  expect_error({
    vecchia_parent_structure(
      data.frame(x = 1 : 5, z = 1 : 5),
      geowarp_model(
        variable = 'value',
        horizontal_coordinates = 'x',
        vertical_coordinate = 'z',
        horizontal_domains = c(0, 10),
        vertical_domain = c(0, 10),
        deviation_model = geowarp_white_deviation_model()
      ),
      n_parents = 5
    )
  })
})

test_that('vecchia_parent_structure with observed only', {
  observed_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  n_parents <- 5

  check_parent_structure <- function(parent_structure) {
    expect_equal(length(parent_structure$observed_ordering), nrow(observed_df))
    expect_equal(nrow(parent_structure$observed_parents), nrow(observed_df))
    expect_equal(ncol(parent_structure$observed_parents), n_parents + 1)
  }
  for (observed_type in c('vertically_dense', 'regular')) {
    parent_structure <- vecchia_parent_structure(
      observed_df,
      full_model,
      n_parents = n_parents,
      observed_type = observed_type
    )
    check_parent_structure(parent_structure)
  }
  parent_structure <- vecchia_parent_structure(
    observed_df,
    vertical_only_model,
    n_parents = n_parents
  )
  check_parent_structure(parent_structure)
})

test_that('vecchia_parent_structure with observed and predicted', {
  observed_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  predicted_df <- data.frame(
    x = rep(0.1 + (1 : 5), each = 5),
    z = rep(0.1 + (1 : 5), 5)
  )
  n_parents <- 10

  check_parent_structure <- function(parent_structure, extra_obs = 0) {
    expect_equal(length(parent_structure$observed_ordering), nrow(observed_df))
    expect_equal(nrow(parent_structure$observed_parents), nrow(observed_df))
    expect_equal(ncol(parent_structure$observed_parents), n_parents + 1)
    expect_equal(length(parent_structure$prediction_ordering), nrow(predicted_df) + extra_obs)
    expect_equal(nrow(parent_structure$prediction_within_parents), nrow(predicted_df) + extra_obs)
    expect_equal(ncol(parent_structure$prediction_within_parents), 6)
    expect_equal(nrow(parent_structure$prediction_between_parents), nrow(predicted_df) + extra_obs)
    expect_equal(ncol(parent_structure$prediction_between_parents), 5)
  }
  for (observed_type in c('vertically_dense', 'regular')) {
    for (predicted_type in c('vertically_dense', 'regular')) {
      parent_structure <- vecchia_parent_structure(
        observed_df,
        full_model,
        n_parents = n_parents,
        observed_type = observed_type,
        prediction_df = predicted_df,
        prediction_type = predicted_type
      )
      check_parent_structure(parent_structure)
    }
  }
  parent_structure <- vecchia_parent_structure(
    observed_df,
    vertical_only_model,
    n_parents = n_parents,
    observed_type = observed_type,
    prediction_df = predicted_df,
    prediction_type = predicted_type
  )
  check_parent_structure(parent_structure)
  parent_structure <- vecchia_parent_structure(
    observed_df,
    vertical_only_model,
    n_parents = n_parents,
    observed_type = observed_type,
    prediction_df = rbind(
      # Including these to target a conditional where some of the prediction
      # locations are the same as the observed locations
      observed_df,
      predicted_df
    ),
    prediction_type = predicted_type
  )
  check_parent_structure(parent_structure, extra_obs = nrow(observed_df))
})

test_that('vecchia_parent_structure with regular scaling works', {
  observed_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  predicted_df <- data.frame(
    x = rep(0.1 + (1 : 5), each = 5),
    z = rep(0.1 + (1 : 5), 5)
  )
  n_parents <- 10

  parent_structure <- vecchia_parent_structure(
    observed_df,
    full_model,
    n_parents = n_parents,
    observed_type = 'regular',
    observed_options = list(scaling = c(0.5, 0.5)),
    prediction_df = predicted_df
  )
  expect_equal(ncol(parent_structure$observed_parents), n_parents + 1)
})

test_that('vecchia_parent_structure with too few neighbours for vertically dense works', {
  observed_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  predicted_df <- data.frame(
    x = rep(0.1 + (1 : 5), each = 5),
    z = rep(0.1 + (1 : 5), 5)
  )
  n_parents <- 4

  expect_warning(
    vecchia_parent_structure(
      observed_df,
      full_model,
      n_parents = n_parents,
      observed_type = 'vertically_dense',
      prediction_df = predicted_df
    ),
    'Too few neighbours compared to horizontal locations; Vecchia approximation may be poor'
  )
})

test_that('vecchia_parent_structure works with partial options', {
  observed_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  predicted_df <- data.frame(
    x = rep(0.1 + (1 : 5), each = 5),
    z = rep(0.1 + (1 : 5), 5)
  )
  n_parents <- 10

  parent_structure <- vecchia_parent_structure(
    observed_df,
    full_model,
    n_parents = n_parents,
    observed_options = list(),
    prediction_df = predicted_df,
    prediction_within_options = list(),
    prediction_between_options = list()
  )
  expect_equal(ncol(parent_structure$observed_parents), n_parents + 1)
})

## vecchia_grouping

test_that('vecchia_grouping has valid output', {
  observed_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  parent_structure <- vecchia_parent_structure(
    observed_df,
    full_model,
    n_parents = 10
  )
  grouping <- vecchia_grouping(parent_structure)

  expect_equal(grouping$N_indices, length(grouping$block_indices))
  expect_equal(grouping$N_blocks, length(grouping$block_last_index))
  expect_equal(grouping$N_blocks, length(grouping$block_N_responses))
})

## plot_parent_structure

test_that('plot_parent_structure works with observed only', {
  observed_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  parent_structure <- vecchia_parent_structure(
    observed_df,
    full_model,
    n_parents = 10
  )
  result <- plot_parent_structure(observed_df, full_model, parent_structure)
  expect_true(is.ggplot(result))
})

test_that('plot_parent_structure works with observed and predicted', {
  base_df <- data.frame(
    x = rep(1 : 5, each = 5),
    z = rep(1 : 5, 5)
  )
  parent_structure <- vecchia_parent_structure(
    base_df,
    full_model,
    n_parents = 10,
    prediction_df = base_df
  )
  result <- plot_parent_structure(
    base_df,
    full_model,
    parent_structure,
    base_df
  )
  expect_true(is.ggplot(result))
})
