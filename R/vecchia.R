#' @export
get_parent_structure <- function(
  df,
  model,
  n_parents = 10,
  scaling = 'auto',
  ordering = GpGp::order_maxmin
) {
  if (is.data.frame(df)) {
    df <- list(df)
  }
  if (scaling[1] == 'auto') {
    log_debug('Computing automatic Vecchia dimensional scaling')
    scaling <- auto_vecchia_scaling(
      df[[1]],
      model,
      n_parents,
      ordering = ordering
    )
  }
  ordering <- match.fun(ordering)
  x_scaled <- lapply(df, function(df_i) {
    x <- as.matrix(
      df_i[, c(model$horizontal_coordinates, model$vertical_coordinate)]
    )
    t(t(x) * scaling)
  })
  log_debug('Finding Vecchia ordering')
  x_ordering <- lapply(x_scaled, ordering)
  log_debug('Finding Vecchia parents')
  parents <- GpGp::find_ordered_nn(
    do.call(rbind, lapply(seq_along(x_scaled), function(i) {
      x_scaled[[i]][x_ordering[[i]], ]
    })),
    n_parents
  )
  list(
    ordering = if (length(x_ordering) == 1) {
      x_ordering[[1]]
    } else {
      x_ordering
    },
    parents = parents
  )
}

#' @export
auto_vecchia_scaling <- function(
  df,
  model,
  n_parents,
  proportion_from_neighbours_target = 0.5,
  tolerance = 0.05,
  max_attempts = 50,
  interval = c(1e-6, 1e5),
  ...
) {
  as_scaling <- function(x) c(rep(1, length(model$horizontal_coordinates)), x)
  bounds <- interval
  previous_proportion_from_neighbours <- 2
  for (i in seq_len(max_attempts)) {
    log_trace('Step {i}: bounds = {bounds[1]}, {bounds[2]}; proportion = {previous_proportion_from_neighbours}')
    middle <- mean(bounds)
    parent_structure_i <- get_parent_structure(
      df,
      model,
      n_parents,
      scaling = as_scaling(middle),
      ...
    )
    n_from_neighbours_i <- parents_n_from_neighbours(df, model, parent_structure_i)
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
  parent_structure
) {
  U_parts <- lapply(seq_len(nrow(parent_structure$parents)), function(k) {
    block_size_k <- min(k, ncol(parent_structure$parents))
    indices_k <- rev(parent_structure$parents[k, 1 : block_size_k])

    Sigma_k <- .covariance_matrix_internal(
      x[indices_k, , drop = FALSE],
      X_deviation_fixed[indices_k, , drop = FALSE],
      X_deviation_random[indices_k, , drop = FALSE],
      fit
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
