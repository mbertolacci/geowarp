#' The 'geowarp' package.
#'
#' @description A package to fit GeoWarp, a model that uses warped spatial
#' processes for inferring soil properties.
#'
#' @docType package
#' @name geowarp-package
#' @aliases geowarp
#' @useDynLib geowarp, .registration = TRUE
#' @import Matrix
#' @import methods
#' @import Rcpp
#' @import dplyr
#' @import ggplot2
#' @importFrom logger log_debug
#' @importFrom logger log_trace
#' @importFrom withr defer
#' @importFrom rstan sampling
#'
NULL

.onLoad <- function(libname, pkgname) {
  base_options <- options()
  geowarp_options <- list(
    geowarp.threads = 1L
  )
  to_set <- !(names(geowarp_options) %in% names(base_options))
  if (any(to_set)) options(geowarp_options[to_set])

  invisible()
}
