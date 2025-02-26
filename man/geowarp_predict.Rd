% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.R
\name{geowarp_predict}
\alias{geowarp_predict}
\alias{predict.geowarp_fit}
\title{Predict from a Fitted GeoWarp Model}
\usage{
geowarp_predict(
  fit,
  prediction_df = fit$observed_df,
  observed_df = fit$observed_df,
  parameters = fit$parameters,
  model = fit$model,
  include_mean = TRUE,
  include_precision_V = FALSE,
  include_samples = FALSE,
  nugget = TRUE,
  n_samples = 500,
  vecchia = "auto",
  n_parents = 50,
  parent_structure = vecchia_parent_structure(observed_df, model, n_parents,
    prediction_df = prediction_df),
  threads = getOption("geowarp.threads"),
  ...
)

\method{predict}{geowarp_fit}(object, ...)
}
\arguments{
\item{fit}{GeoWarp model fit produced using \code{\link{geowarp_optimise}}.}

\item{prediction_df}{Data frame for locations where predictions are to be
made. Defaults to the observed data frame in the model fit.}

\item{observed_df}{Data frame containing the observed data. Defaults to the}

\item{parameters}{GeoWarp parameters (default uses the parameters from the
fit).}

\item{model}{GeoWarp model object (default uses the model from the fit).}

\item{include_mean}{Logical, indicating whether to include the mean in the
output.}

\item{include_precision_V}{Logical, indicating whether to include the V
matrix in the output (explained below).}

\item{include_samples}{Logical, indicating whether to include the predictive
samples in the output.}

\item{nugget}{Logical, indicating whether to include a nugget effect in the
prediction.}

\item{n_samples}{Number of predictive samples to generate.}

\item{vecchia}{Whether to use Vecchia's approximation (\code{'auto'}, \code{TRUE},
\code{FALSE}). If \code{'auto'}, the Vecchia approximation is used when the number of
observed and prediction locations together exceed 1000.}

\item{n_parents}{Number of parents to consider for Vecchia's approximation.}

\item{parent_structure}{Parent structure for Vecchia's approximation,
produced using \code{\link{vecchia_parent_structure}}.}

\item{threads}{Number of threads to use for parallelisation, taken from
\code{getOption('geowarp.threads')} by default, which is itself set to 1 by
default. The special value -1 picks a number of threads based on the number
of cores in the system). Only applies if using the Vecchia approximation.}

\item{...}{Ignored.}
}
\value{
A list containing the elements of the predictive distribution
(depending on the options):
\itemize{
\item \code{mean}: The mean of the predictive distribution.
\item \code{precision_V}: A (possibly sparse) matrix V such that VV' is the
reordered precision matrix of the predictive distribution.
\item \code{ordering}: A vector of integers indicating the ordering of the
entries in the matrix V; see Details below.
\item \code{samples}: A matrix of predictive samples.
}
}
\description{
This function generates predictions using a fitted GeoWarp model. The
prediction can use the Vecchia approximation when the number of observed and
predicted data are large. The prediction distribution is multivariate
Gaussian, so it is specified by its mean and covariance
}
\section{Functions}{
\itemize{
\item \code{predict(geowarp_fit)}: Predict using a fit object

}}
\section{Details}{

If \code{prediction} contains the output of this function, the precision matrix
of the prediction can be recovered by calculating
\code{tcrossprod(prediction$precision_V)[prediction$ordering, prediction$ordering]}.
However, this may be of limited use because, while the precision matrix may
be sparse, the covariance may not be. To calculate variances and other
quantities it is better to draw samples from the prediction distribution.
This function can do that, or you can do it manually as shown in the example
section.
}

\examples{
# Assuming `fit` is a fitted GeoWarp model and `new_data` is a new data frame
# for prediction:
prediction_dist <- predict(
  fit,
  prediction_df = new_data,
  include_precision_V = TRUE,
  include_samples = TRUE
)
# Samples are now in `prediction_dist$samples`, but you can also generate
# them manually
y_tilde <- as.matrix(solve(
  t(prediction_dist$precision_V),
  z
))
manual_samples <- prediction_dist$mean + y_tilde[prediction_dist$ordering, ]

}
\seealso{
\code{\link{geowarp_optimise}}
\code{\link{vecchia_parent_structure}}
}
