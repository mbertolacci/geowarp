#' Calculate a Parent Structure for the Vecchia Approximation
#'
#' This function calculates the parent structure for both observed and optional
#' predicted locations based on the Vecchia approximation. It's used to
#' structure the covariance or conditional dependencies in GeoWarp for both
#' inference and prediction. All observed locations are placed earlier in
#' the ordering than the predicted locations.
#'
#' The chosen parent structure depends on the pattern of the observed and
#' prediction locations, which is user specified. The observed locations are
#' conditioned on each other, while the predicted locations condition on both
#' the observed and the other predicted locations to a degree specified by the
#' user.
#'
#' @param observed_df Data frame containing observed locations.
#' @param model GeoWarp model objects.
#' @param n_parents Number of parents to use. Default is 20.
#' @param observed_ordering Ordering of the observed locations. Defaults to
#' a random reordering, except for vertical only models.
#' @param observed_type Pattern of the observed locations. One of
#' 'vertically_dense' (the default, useful for CPT data), or 'regular' (useful
#' for data on a grid or that is uniformly distributed).
#' @param observed_options List of additional options for observed structure.
#' See details below.
#' @param prediction_df Optional data frame containing prediction locations.
#' @param prediction_type Pattern of the predicted locations; see
#' \code{observed_type}. Defaults to regular, for a grid.
#' @param prediction_ordering Ordering of the prediction locations; see
#' `observed_ordering`.
#' @param prediction_within_options Options for the parent structure within
#' the prediction locations.
#' @param prediction_between_options Options for the parent structure between
#' the prediction and the observation locations.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{observed_ordering}: A vector of indices specifying the ordering
#' used for observed locations.
#' \item \code{observed_parents}: The parent structure for observed locations.
#' Given as a matrix. The first column contains the index of the (reordered)
#' location, and the subsequent columns contain the indices of the parents (some
#' of which can be `NA`).
#' \item \code{prediction_ordering}: The ordering used for prediction locations
#' (if provided).
#' \item \code{prediction_within_parents}: Matrix of parents for within the
#' prediction locations (if provided).
#' \item \code{prediction_between_parents}: Matrix of parents for between the
#' predictions and observations.
#' }
#'
#' @section Location pattern and chosen parent structure:
#' The observed data condition on themselves. For
#' \code{observed_type = 'vertically_dense'}, the scheme is to split the
#' parents, by default 50-50, between the nearest neighbours of the points among
#' its predecessors in the ordering, and points at other horizontal locations.
#' This split can be changed by setting `observed_options$n_parents` and
# `observed_options$n_horizontal_parents`.
#' For \code{observed_type = 'regular'}, the parents are just the nearest
#' neighbours among the predecessor points. Then `observed_options` must contain
#' the number of parents in `n_parents` and, optionally, a vector named
#' `scaling` by which each coordinate will be rescaled; this can be used if the
#' grid is denser in one direction than another to ensure the neighbourhood
#' structure is regular too.
#'
#' The prediction locations condition on both the observed locations and the
#' other prediction locations. The options for the scheme for conditioning
#' within the predictions is the same as for the observed locations, and the
#' behaviour of this scheme is controlled through `prediction_within_options`,
#' which has the same options as `observed_options`.
#'
#' For the between-parents, the scheme depends on `observed_type`: if the
#' observations are vertically dense, the observation parents are spread among
#' the different observed horizontal locations; otherwise, the nearest
#' neighbours among the parents are used. In both cases, the number of between
#' parents is controlled by `prediction_between_options$n_parents`. For
#' `observed_type = 'regular'`, the `scaling` option can be used for the
#' same effect as for the observed locations.
#' @export
vecchia_parent_structure <- function(
  observed_df,
  model,
  n_parents = 20,
  observed_ordering,
  observed_type = c('vertically_dense', 'regular'),
  observed_options = list(n_parents = n_parents),
  prediction_df,
  prediction_type = c('regular', 'vertically_dense'),
  prediction_ordering,
  prediction_within_options = list(n_parents = floor(n_parents / 2)),
  prediction_between_options = list(n_parents = ceiling(n_parents / 2))
) {
  observed_type <- match.arg(observed_type)
  prediction_type <- match.arg(prediction_type)

  if (model$deviation_model$name == 'white') {
    stop('Vecchia not supported (or needed) for white covariance')
  }

  if (is.null(observed_options$n_parents)) {
    observed_options$n_parents <- n_parents
  }

  # Within observed structure
  observed_structure <- .parent_structure_within(
    observed_df,
    model,
    observed_ordering,
    observed_type,
    observed_options
  )
  output <- list(
    observed_ordering = observed_structure$ordering,
    observed_parents = observed_structure$parents
  )

  if (!missing(prediction_df)) {
    if (is.null(prediction_within_options$n_parents)) {
      prediction_within_options$n_parents <- floor(n_parents / 2)
    }
    if (is.null(prediction_between_options$n_parents)) {
      prediction_between_options$n_parents <- ceiling(n_parents / 2)
    }

    # Prediction within structure
    prediction_within_structure <- .parent_structure_within(
      prediction_df,
      model,
      prediction_ordering,
      prediction_type,
      prediction_within_options
    )

    # Prediction between structure
    prediction_between_method <- if (model$deviation_model$name == 'vertical_only') {
      .between_parents_vertical_only
    } else if (observed_type == 'vertically_dense') {
      .between_parents_vertically_dense
    } else if (observed_type == 'regular') {
      .between_parents_regular
    }

    prediction_between_parents <- prediction_between_method(
      observed_df,
      prediction_df,
      model,
      observed_structure$ordering,
      prediction_within_structure$ordering,
      prediction_between_options
    )

    output$prediction_ordering <- prediction_within_structure$ordering
    output$prediction_within_parents <- prediction_within_structure$parents
    output$prediction_between_parents <- prediction_between_parents
  }

  output
}

#' Calculate Groups of Locations
#'
#' Takes a parent structure and automatically partition the array into groups
#' that share many parents. In each group, some observations are conditioned on,
#' and some are the responses. No response ever belongs to more than one group.
#' This is used in GeoWarp's model fitting to calculate the model likelihood.
#' Under this hood this uses \code{\link[GpGp]{group_obs}} function in the GpGp
#' package.
#'
#' @param parent_structure A parent structure list as returned by
#' \code{\link{vecchia_parent_structure}}.
#' @param exponent Criterion for grouping, passed to `GpGp::group_obs`. Should
#' be either 2 or 3.
#' @return A list with the following entries:
#' \itemize{
#' \item \code{N_indices}: The length of `all_inds`.
#' \item \code{N_blocks}: The number of groups/blocks.
#' \item \code{block_indices}: A vector of which indices belong to which group.
#' Contains the indices in block 1 first, followed by those in block 2, and so
#' on. The indices in each block are ordered with the responses at the end.
#' \item \code{block_last_index}: For each block, the last index in
#' `block_indices` that contains the block's data indices.
#' \item \code{block_N_responses}: For each block, the number of responses in
#' the block.
#' }
#' @seealso
#' \code{\link[GpGp]{group_obs}}
#' @export
vecchia_grouping <- function(
  parent_structure,
  exponent = 3
) {
  vecchia_groups <- GpGp::group_obs(
    parent_structure$observed_parents,
    exponent = exponent
  )
  current_block_index <- 1L
  current_response_index <- 1L
  for (i in seq_along(vecchia_groups$last_ind_of_block)) {
    block_i <- current_block_index : vecchia_groups$last_ind_of_block[i]
    block_indices <- vecchia_groups$all_inds[block_i]
    response_indices <- vecchia_groups$local_resp_inds[
      current_response_index : vecchia_groups$last_resp_of_block[i]
    ]
    vecchia_groups$all_inds[block_i] <- c(
      block_indices[-response_indices],
      block_indices[response_indices]
    )
    current_block_index <- vecchia_groups$last_ind_of_block[i] + 1L
    current_response_index <- vecchia_groups$last_resp_of_block[i] + 1L
  }
  list(
    N_indices = length(vecchia_groups$all_inds),
    N_blocks = length(vecchia_groups$last_ind_of_block),
    block_indices = vecchia_groups$all_inds,
    block_last_index = vecchia_groups$last_ind_of_block,
    block_N_responses = c(
      vecchia_groups$last_resp_of_block[1],
      diff(vecchia_groups$last_resp_of_block)
    )
  )
}

#' Plot a Vecchia Parent Structure
#'
#' This function visualizes a Vecchia parent structure of data in a 2D plot.
#' Each panel in the plot is either the target, one of its parents, or other.
#'
#' @param observed_df A data frame containing the observed data.
#' @param model GeoWarp model object.
#' @param structure A list containing parent structures and orderings produced
#' by \code{\link{vecchia_parent_structure}}.
#' @param prediction_df A data frame containing the prediction data (optional).
#' @param observed_indices Indices of observed data points to plot (optional;
#' defaults to three chosen at random).
#' @param prediction_indices Indices of prediction data points to plot
#' (optional; defaults to three chosen at random).
#' @param left_coordinate Coordinate to use for the plot x-axis (default is 1).
#' @param right_coordinate Coordinate to use for the plot y-axis (default is 2).
#' @return A ggplot2 plot.
#' @export
plot_parent_structure <- function(
  observed_df,
  model,
  structure,
  prediction_df,
  observed_indices,
  prediction_indices,
  left_coordinate = 1,
  right_coordinate = 2
) {
  .get_plot_df <- function(x, indices, parents, type, is_within, prefix) {
    n <- nrow(x)
    lapply(indices, function(i) {
      tibble(
        index = sprintf('%s%d', prefix, i),
        type = type,
        left = x[, left_coordinate],
        right = x[, right_coordinate],
        status = factor(case_when(
          (1 : n) == i & is_within ~ 'target',
          (1 : n) %in% parents[i, ] ~ 'parent',
          TRUE ~ 'other'
        ), c('other', 'parent', 'target'))
      )
    }) %>%
      bind_rows()
  }

  if (missing(observed_indices)) {
    observed_indices <- sample.int(nrow(observed_df), 3L)
  }

  x_observed <- model_coordinates(observed_df, model)[
    structure$observed_ordering,
  ]

  plot_observed_df <- .get_plot_df(
    x_observed,
    observed_indices,
    structure$observed_parents,
    'observed',
    TRUE,
    'observed'
  )

  plot_prediction_df <- NULL
  if (!missing(prediction_df)) {
    if (missing(prediction_indices)) {
      prediction_indices <- sample.int(nrow(prediction_df), 3)
    }

    x_prediction <- model_coordinates(prediction_df, model)[
      structure$prediction_ordering,
    ]
    plot_prediction_df <- bind_rows(
      .get_plot_df(
        x_observed,
        prediction_indices,
        structure$prediction_between_parents,
        'observed',
        FALSE,
        'prediction'
      ),
      .get_plot_df(
        x_prediction,
        prediction_indices,
        structure$prediction_within_parents,
        'prediction',
        TRUE,
        'prediction'
      )
    )
  }

  bind_rows(plot_observed_df, plot_prediction_df) %>%
    arrange(status) %>%
    ggplot(aes(left, right, colour = status, shape = type)) +
      geom_point() +
      scale_colour_manual(values = c(
        'target'= 'red',
        'parent' = 'blue',
        'other' = '#aaaaaa'
      )) +
      facet_wrap(~ index)
}

.parent_structure_within <- function(
  df,
  model,
  ordering,
  type,
  options
) {
  is_vertical_only <- model$deviation_model$name == 'vertical_only'

  if (is_vertical_only) {
    # Ignore any provided
    ordering <- .vertical_only_ordering(df, model)
  } else if (missing(ordering)) {
    ordering <- sample.int(nrow(df))
  }

  # Within observed structure
  parents <- if (is_vertical_only) {
    .parents_vertical_only(
      df,
      model,
      ordering,
      options
    )
  } else if (type == 'vertically_dense') {
    .parents_vertically_dense(
      df,
      model,
      ordering,
      options
    )
  } else if (type == 'regular') {
    .parents_regular(
      df,
      model,
      ordering,
      options
    )
  }

  list(
    ordering = ordering,
    parents = parents
  )
}

.vertical_only_ordering <- function(df, model) {
  do.call(order, lapply(
    c(model$horizontal_coordinates, model$vertical_coordinate),
    function(name) df[[name]]
  ))
}

.parents_vertical_only <- function(
  df,
  model,
  ordering,
  options
) {
  n_parents <- options$n_parents

  x <- model_coordinates(df, model)[ordering, ]

  vertical_index <- ncol(x)
  horizontal_indices <- seq_len(ncol(x) - 1)
  horizontal_groups <- attr(.unique_rows(
    x[, horizontal_indices, drop = FALSE]
  ), 'indices')

  parent_parts <- lapply(split(
    seq_len(nrow(x)),
    horizontal_groups
  ), function(indices) {
    output <- parent_knn(
      x[indices, vertical_index],
      n_parents
    )
    # Use the canonical indices derived above in the parent object
    output[] <- indices[output]
    output
  })

  parents <- do.call(rbind, parent_parts)
  parents[order(parents[, 1]), ]
}

.parents_vertically_dense <- function(
  df,
  model,
  ordering,
  options
) {
  if (is.null(options$n_nearest_parents)) {
    options$n_nearest_parents <- floor(options$n_parents / 2)
  }
  if (is.null(options$n_horizontal_parents)) {
    options$n_horizontal_parents <- options$n_parents - options$n_nearest_parents
  }

  x <- model_coordinates(df, model)[ordering, ]

  parents_nearest <- parent_knn(x, options$n_nearest_parents)

  horizontal_groups <- attr(.unique_rows(
    x[, 1 : (ncol(x) - 1), drop = FALSE]
  ), 'indices')
  parents_horizontal <- parent_knn_across_groups(
    x[, ncol(x)],
    horizontal_groups,
    options$n_horizontal_parents
  )
  # Drop redundant first column
  parents_horizontal <- parents_horizontal[ , 2 : ncol(parents_horizontal)]

  .deduplicate_parents(cbind(parents_nearest, parents_horizontal))
}

.parents_regular <- function(
  df,
  model,
  ordering,
  options
) {
  x <- model_coordinates(df, model)[ordering, ]
  if (!is.null(options$scaling)) {
    x <- t(options$scaling * t(x))
  }
  parent_knn(x, options$n_parents)
}

.between_parents_vertical_only <- function(
  observed_df,
  prediction_df,
  model,
  observed_ordering,
  prediction_ordering,
  options
) {
  observed_x <- model_coordinates(observed_df, model)[observed_ordering, ]
  prediction_x <- model_coordinates(prediction_df, model)[prediction_ordering, ]

  vertical_index <- ncol(observed_x)
  horizontal_indices <- seq_len(ncol(observed_x) - 1)
  horizontal_groups <- attr(.unique_rows(
    rbind(
      observed_x[, horizontal_indices, drop = FALSE],
      prediction_x[, horizontal_indices, drop = FALSE]
    )
  ), 'indices')
  first_prediction_index <- nrow(observed_x) + 1

  all_parents <- lapply(split(
    seq_len(nrow(observed_x) + nrow(prediction_x)),
    horizontal_groups
  ), function(indices) {
    observed_indices <- indices[indices < first_prediction_index]
    prediction_indices <- (
      indices[indices >= first_prediction_index]
      - first_prediction_index
      + 1L
    )

    if (length(prediction_indices) == 0) {
      return(NULL)
    }

    output_i <- matrix(
      NA_integer_,
      nrow = length(prediction_indices),
      ncol = options$n_parents + 1
    )
    output_i[, 1] <- prediction_indices
    if (length(observed_indices) > 0) {
      output_i[, 2 : (options$n_parents + 1)] <- .knn(
        observed_x[observed_indices, vertical_index, drop = FALSE],
        prediction_x[prediction_indices, vertical_index, drop = FALSE],
        options$n_parents
      )
    }
    output_i
  })
  parents <- do.call(rbind, all_parents)
  # Order by the first column, which holds the index
  parents <- parents[order(parents[, 1]), ]

  parents[, 2 : ncol(parents), drop = FALSE]
}

# Scheme: grab the closest neighbours from each horizontal location in the
# observations
.between_parents_vertically_dense <- function(
  observed_df,
  prediction_df,
  model,
  observed_ordering,
  prediction_ordering,
  options
) {
  observed_x <- model_coordinates(observed_df, model)[observed_ordering, ]
  prediction_x <- model_coordinates(prediction_df, model)[prediction_ordering, ]

  vertical_index <- ncol(observed_x)
  horizontal_indices <- seq_len(ncol(observed_x) - 1)
  horizontal_groups <- .unique_rows(
    observed_x[, horizontal_indices, drop = FALSE]
  )

  n_groups <- nrow(horizontal_groups)
  n_parents_per_group <- c(
    rep(floor(options$n_parents / n_groups) + 1, options$n_parents %% n_groups),
    rep(floor(options$n_parents / n_groups), n_groups - options$n_parents %% n_groups)
  )

  if (any(n_parents_per_group == 0)) {
    warning('Too few neighbours compared to horizontal locations; Vecchia approximation may be poor')
  }

  output <- matrix(
    NA_integer_,
    nrow = nrow(prediction_x),
    ncol = options$n_parents
  )
  for (group_i in seq_len(n_groups)) {
    if (n_parents_per_group[group_i] == 0) next

    start_index <- if (group_i == 1) {
      1
    } else {
      cumsum(n_parents_per_group)[group_i - 1] + 1
    }

    output_indices <- start_index : (start_index + n_parents_per_group[group_i] - 1)
    observed_indices_i <- which(attr(horizontal_groups, 'indices') == group_i)
    parents_i <- .knn(
      observed_x[observed_indices_i, vertical_index, drop = FALSE],
      prediction_x[, vertical_index, drop = FALSE],
      n_parents_per_group[group_i]
    )
    parents_i[] <- observed_indices_i[parents_i]

    output[, output_indices] <- parents_i
  }

  output
}

.between_parents_regular <- function(
  observed_df,
  prediction_df,
  model,
  observed_ordering,
  prediction_ordering,
  options
) {
  x1 <- model_coordinates(observed_df, model)[observed_ordering, ]
  x2 <- model_coordinates(prediction_df, model)[prediction_ordering, ]
  if (!is.null(options$scaling)) {
    x1 <- t(options$scaling * t(x1))
    x2 <- t(options$scaling * t(x2))
  }
  .knn(x1, x2, options$n_parents)
}

.create_vecchia_U <- function(
  x,
  X_deviation_fixed,
  X_deviation_random,
  fit,
  parents,
  nugget = rep(TRUE, nrow(x)),
  model = fit$model,
  parameters = fit$parameters,
  threads
) {
  x_warped <- warped_coordinates(x = x, model = model, parameters = parameters)
  if (model$deviation_model$name == 'vertical_only') {
    # NOTE(mgnb): The neighbourhood structure itself will never pick a parent
    # from a different vertical location, so this is safe
    x_warped <- x_warped[, ncol(x_warped), drop = FALSE]
  }

  sigma_deviation <- exp(
    0.5 * as.vector(
      X_deviation_fixed %*% parameters$eta_deviation
      + X_deviation_random %*% parameters$zeta_deviation
    )
  )
  parents_na0 <- parents[, ncol(parents) : 1]
  parents_na0[is.na(parents_na0)] <- 0
  parents_na0 <- rbind(
    rep(0, ncol(parents)),
    parents_na0
  )

  i_entries_raw <- rep(0 : (nrow(parents)), each = ncol(parents))
  j_entries_raw <- as.vector(t(parents_na0))
  x_entries_raw <- vecchia_U_x_parts(
    parents_na0,
    rbind(
      rep(0, ncol(x_warped)),
      x_warped
    ),
    model$deviation_model$covariance_function,
    c(0, sigma_deviation),
    c(TRUE, nugget) * parameters$sigma_squared_nugget,
    threads > 1 || threads == -1
  )

  i_entries <- i_entries_raw[t(parents_na0) != 0]
  j_entries <- j_entries_raw[t(parents_na0) != 0]
  x_entries <- x_entries_raw[t(parents_na0) != 0]

  t(Matrix::sparseMatrix(
    i = i_entries,
    j = j_entries,
    x = x_entries,
    triangular = TRUE,
    dims = c(nrow(x), nrow(x))
  ))
}

.knn <- function(x, y, k) {
  if (k < nrow(x)) {
    FNN::get.knnx(x, y, k)$nn.index
  } else {
    output <- matrix(NA, nrow = nrow(x), ncol = k)
    for (i in seq_len(nrow(x))) {
      output[, i] <- i
    }
    output
  }
}
