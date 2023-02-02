#' @export
pcpt_cross_validation_fit <- function(
  df,
  grouping_variable,
  model,
  method = 'pcpt_optimise',
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
      pbmcapply::pbmclapply(..., ignore.interactive = TRUE, mc.cores = grouping_cores, mc.preschedule = FALSE)
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
  any_error <- FALSE
  for (result in output) {
    if (is(result, 'try-error')) {
      print(result)
      any_error <- TRUE
    }
  }
  if (any_error) {
    stop('cross validation fits failed')
  }
  output
}
