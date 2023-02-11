#' @export
vecchia_parent_structure <- function(
  df,
  model,
  n_parents = 10,
  scaling = 'auto',
  ordering = GpGp::order_maxmin
) {
  if (is.data.frame(df)) {
    df <- list(df)
  }
  if (model$deviation_model$name == 'vertical_only') {
    .vecchia_parent_structure_vertical_only(df, model, n_parents)
  } else if (model$deviation_model$name == 'white') {
    stop('Vecchia not supported (or needed) for white covariance')
  } else {
    .vecchia_parent_structure_full(df, model, n_parents, ordering, scaling)
  }
}

.vecchia_parent_structure_vertical_only <- function(df, model, n_parents) {
  orderings <- lapply(df, function(df_i) {
    do.call(order, lapply(
      c(model$horizontal_coordinates, model$vertical_coordinate),
      function(name) df_i[[name]]
    ))
  })

  coordinates_joined <- lapply(seq_along(df), function(i) {
    output <- df[[i]] %>%
      ungroup() %>%
      select(model$horizontal_coordinates, model$vertical_coordinate)
    output[orderings[[i]], ]
  }) %>%
    bind_rows() %>%
    mutate(index = seq_len(n()))

  parent_parts <- coordinates_joined %>%
    group_by(across(model$horizontal_coordinates)) %>%
    group_map(~ {
      output <- GpGp::find_ordered_nn(.x[[model$vertical_coordinate]], n_parents)
      if (ncol(output) < n_parents + 1) {
        output <- cbind(
          output,
          matrix(NA, nrow = nrow(output), ncol = n_parents + 1 - ncol(output))
        )
      }
      # Use the canonical indices derived above in the parent object
      output[] <- .x$index[output]
      output
    })
  parents <- do.call(rbind, parent_parts)
  parents <- parents[order(parents[, 1]), ]

  list(
    ordering = if (length(orderings) == 1) {
      orderings[[1]]
    } else {
      orderings
    },
    parents = parents
  )
}

.vecchia_parent_structure_full <- function(df, model, n_parents, ordering, scaling) {
  if (!is.list(scaling) && scaling == 'auto') {
    scaling <- rep(list('auto'), length(df))
  }
  if (!is.list(scaling)) {
    scaling <- list(scaling)
  }
  stopifnot(length(scaling) == length(df))
  scaling <- lapply(seq_along(scaling), function(i) {
    scaling_i <- scaling[[i]]
    if (length(scaling_i) > 1 && scaling_i != 'auto') return(scaling[[i]])

    log_debug('Computing automatic Vecchia dimensional scaling')
    auto_vecchia_scaling(
      df[[i]],
      model,
      n_parents,
      parent_df = head(df, i - 1),
      ordering = ordering
    )
  })
  ordering <- match.fun(ordering)
  x_parts <- lapply(seq_along(df), function(i) {
    as.matrix(
      df[[i]][, c(model$horizontal_coordinates, model$vertical_coordinate)]
    )
  })

  orderings <- lapply(seq_along(x_parts), function(i) {
    ordering(t(t(x_parts[[i]]) * scaling[[i]]))
  })
  x_ordered <- do.call(rbind, lapply(seq_along(x_parts), function(i) {
    x_parts[[i]][orderings[[i]], ]
  }))

  part_sizes <- sapply(df, nrow)
  parents <- do.call(rbind, lapply(seq_along(df), function(i) {
    max_index <- sum(part_sizes[seq_len(i)])
    x_ordered_upto_scaled <- t(t(x_ordered[seq_len(max_index), ]) * scaling[[i]])
    tail(
      GpGp::find_ordered_nn(
        x_ordered_upto_scaled,
        n_parents
      ),
      part_sizes[i]
    )
  }))

  list(
    ordering = if (length(orderings) == 1) {
      orderings[[1]]
    } else {
      orderings
    },
    parents = parents,
    scaling = scaling
  )
}

#' @export
vecchia_grouping <- function(
  parent_structure,
  exponent = 2
) {
  vecchia_groups <- GpGp::group_obs(
    parent_structure$parents,
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

#' @export
auto_vecchia_scaling <- function(
  df,
  model,
  n_parents,
  parent_df,
  proportion_from_neighbours_target = 0.5,
  tolerance = 0.05,
  max_attempts = 50,
  interval = c(1e-6, 1e5),
  ...
) {
  if (missing(parent_df)) {
    parent_df <- NULL
  }

  if (is.null(parent_df) && nrow(distinct(df[, model$horizontal_coordinates, drop = FALSE])) == 1) {
    return(rep(1, 1 + length(model$horizontal_coordinates)))
  }

  as_scaling <- function(x) c(rep(1, length(model$horizontal_coordinates)), x)
  bounds <- interval
  previous_proportion_from_neighbours <- 2
  for (i in seq_len(max_attempts)) {
    log_trace('Step {i}: bounds = {bounds[1]}, {bounds[2]}; proportion = {previous_proportion_from_neighbours}')
    middle <- mean(bounds)
    parent_structure_i <- vecchia_parent_structure(
      c(parent_df, list(df)),
      model,
      n_parents,
      scaling = rep(
        list(as_scaling(middle)),
        length(parent_df) + 1
      ),
      ...
    )
    n_from_neighbours_i <- parents_n_from_neighbours(
      df,
      model,
      parent_structure_i,
      parent_df
    )
    proportion_from_neighbours_i <- mean(n_from_neighbours_i) / n_parents
    if (abs(proportion_from_neighbours_i - proportion_from_neighbours_target) <= tolerance) {
      return(as_scaling(middle))
    }

    if (proportion_from_neighbours_i < proportion_from_neighbours_target) {
      bounds[1] <- middle
    } else if (proportion_from_neighbours_i > proportion_from_neighbours_target) {
      bounds[2] <- middle
    } else {
      return(as_scaling(middle))
    }
    previous_proportion_from_neighbours <- proportion_from_neighbours_i
  }
  stop('max attempts exceeded')
}

#' @export
parents_n_from_neighbours <- function(
  df,
  model,
  structure,
  parent_df
) {
  if (missing(parent_df)) {
    parent_df <- NULL
  }
  if (length(parent_df) == 0) {
    structure$ordering <- list(structure$ordering)
  }

  all_df <- c(parent_df, list(df))
  x <- do.call(rbind, lapply(seq_along(all_df), function(i) {
    as.matrix(all_df[[i]][, model$horizontal_coordinates])[structure$ordering[[i]], , drop = FALSE]
  }))
  apply(tail(structure$parents, nrow(df)), 1, function(parents_i) {
    x_i <- x[parents_i[1], ]
    output_i <- 0L
    for (k in tail(parents_i, -1)) {
      if (is.na(k)) break
      if (any(x_i != x[k, ])) {
        output_i <- output_i + 1L
      }
    }
    output_i
  })
}

.parents_n_from_parents <- function(
  df,
  parent_df,
  model,
  structure
) {
  x <- as.matrix(df[, model$horizontal_coordinates])[structure$ordering, , drop = FALSE]
  apply(structure$parents, 1, function(parents_i) {
    x_i <- x[parents_i[1], ]
    output_i <- 0L
    for (k in tail(parents_i, -1)) {
      if (is.na(k)) break
      if (any(x_i != x[k, ])) {
        output_i <- output_i + 1L
      }
    }
    output_i
  })
}

#' @export
plot_parent_structure <- function(
  df,
  model,
  structure,
  scaling = 1,
  indices,
  left_coordinate = 1,
  right_coordinate = 2
) {
  if (is.data.frame(df)) {
    df <- list(df)
    structure$ordering <- list(structure$ordering)
  }
  if (missing(indices)) {
    indices <- sample.int(sum(sapply(df, nrow)), 5)
  }

  x <- do.call(rbind, lapply(seq_along(df), function(i) {
    x_i <- as.matrix(
      df[[i]][, c(model$horizontal_coordinate, model$vertical_coordinate)]
    )[structure$ordering[[i]], ]
    t(t(x_i) * scaling)
  }))

  bind_rows(lapply(indices, function(i) {
    n <- nrow(x)
    tibble(
      index = i,
      left = x[, left_coordinate],
      right = x[, right_coordinate],
      status = factor(case_when(
        (1 : n) == i ~ 'target',
        (1 : n) %in% structure$parents[i, ] ~ 'parent',
        TRUE ~ 'other'
      ), c('other', 'parent', 'target'))
    )
  })) %>%
    arrange(status) %>%
    ggplot(aes(left, right, colour = status, shape = status)) +
      geom_point() +
      scale_colour_manual(values = c(
        'target'= 'red',
        'parent' = 'blue',
        'other' = '#aaaaaa'
      )) +
      facet_wrap(~ index)
}

.create_vecchia_U <- function(
  x,
  X_deviation_fixed,
  X_deviation_random,
  fit,
  parent_structure,
  model = fit$model,
  parameters = fit$parameters
) {
  U_parts <- lapply(seq_len(nrow(parent_structure$parents)), function(k) {
    block_size_k <- sum(!is.na(parent_structure$parents[k, ]))
    indices_k <- rev(parent_structure$parents[k, 1 : block_size_k])

    Sigma_k <- .covariance_matrix_internal(
      x[indices_k, , drop = FALSE],
      X_deviation_fixed[indices_k, , drop = FALSE],
      X_deviation_random[indices_k, , drop = FALSE],
      model = model,
      parameters = parameters
    )
    R_k <- chol(Sigma_k)
    list(
      i = rep(k, block_size_k),
      j = indices_k,
      x = backsolve(R_k, c(rep(0, block_size_k - 1), 1))
    )
  })
  t(Matrix::sparseMatrix(
    i = do.call(c, lapply(U_parts, getElement, 'i')),
    j = do.call(c, lapply(U_parts, getElement, 'j')),
    x = do.call(c, lapply(U_parts, getElement, 'x')),
    triangular = TRUE,
    dims = c(nrow(x), nrow(x))
  ))
}
