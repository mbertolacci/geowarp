#' Perform Cross-Validation for GeoWarp
#'
#' This function performs k-fold cross-validation on a specified dataset and
#' GeoWarp model. Each group defined by `grouping_variable` is left out once
#' while fitting the model with the remaining data. The function returns a list
#' of fitted models for each left-out group.
#'
#' @param df Data frame containing the observations.
#' @param grouping_variable Variable that defines the grouping in `df`.
#' @param model A GeoWarp model made using \code{\link{geowarp_model}} that
#' specifies the model structure and priors.
#' @param method The method used for fitting the model, defaults to
#' \code{\link{geowarp_optimise}}.
#' @param show_progress Whether to show progress, defaults to FALSE.
#' @param grouping_cores Number of fits to run in parallel.
#' @param group_output_pattern A \code{\link[glue]{glue}} pattern for naming
#' output files used for logging output (optional).
#' @param prepare_input A function to prepare the input data for each group,
#' defaults to identity function.
#' @param ... Additional arguments to be passed to the method for fitting.
#'
#' @return A list containing the cross-validated fits for each group.
#'
#' @examples
#' # Using geowarp_cross_validation_fit
#' cv_fits <- geowarp_cross_validation_fit(
#'   df = my_data,
#'   grouping_variable = 'group_id',
#'   model = my_model,
#'   show_progress = TRUE
#' )
#'
#' @export
geowarp_cross_validation_fit <- function(
  df,
  grouping_variable,
  model,
  method = 'geowarp_optimise',
  show_progress = FALSE,
  grouping_cores = 1,
  group_output_pattern,
  prepare_input = function(input_df, group) input_df,
  ...
) {
  method <- match.fun(method)

  groups <- df[[grouping_variable]]

  run_jobs <- function(...) {
    if (show_progress) {
      pbmcapply::pbmclapply(..., mc.cores = grouping_cores, mc.preschedule = FALSE)
    } else {
      parallel::mclapply(..., mc.cores = grouping_cores, mc.preschedule = FALSE)
    }
  }

  write_output <- !missing(group_output_pattern)

  output <- run_jobs(
    unique(groups),
    function(group) {
      df_i <- prepare_input(
        df[df[[grouping_variable]] != group, ],
        group
      )
      if (write_output) {
        sink(glue::glue(group_output_pattern))
      }
      output <- list(
        group = group,
        fit = method(
          df = df_i,
          model = model,
          ...
        )
      )
      if (write_output) {
        sink(NULL)
      }
      output
    }
  )

  stopifnot(all(!sapply(output, is, 'try-error')))

  output
}
