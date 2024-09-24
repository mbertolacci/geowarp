.with_tbb_threads <- withr::with_(.set_tbb_threads)
.local_tbb_threads <- withr::local_(.set_tbb_threads)

#' Get and set number of TBB threads
#'
#' Get and set the number of threads used by Intel's Threading Building Blocks
#' (TBB) library, used within Stan for parallelisation.
#'
#' @param n Number of threads
#' @return The number of threads prior to the call
#' @export
set_tbb_threads <- function(n) {
  .set_tbb_threads(n)
}

#' @describeIn set_tbb_threads Get number of TBB threads
#' @export
get_tbb_threads <- function() {
  .get_tbb_threads()
}
