% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RasperGade_alias.R,
%   R/RasperGade_cross_validation.R
\name{crossValidationWithPE}
\alias{crossValidationWithPE}
\alias{LOO_CV_with_PE}
\title{Leave-one-out cross-validation under pulsed evolution model}
\usage{
crossValidationWithPE(
  FMR,
  add.epsilon = TRUE,
  laplace = FALSE,
  numApprox = 1,
  margin = 1e-06,
  numCores = 1,
  asymptotic = 5
)

LOO_CV_with_PE(
  FMR,
  add.epsilon = TRUE,
  laplace = FALSE,
  numApprox = 1,
  margin = 1e-06,
  numCores = 1,
  asymptotic = 5
)
}
\arguments{
\item{FMR}{the returned data structure from function \code{fullMarginalReconstructionWithPE}}

\item{add.epsilon}{logical, if true, time-independent variation is added to the variance in the data frame}

\item{laplace}{logical, if true, Laplace distribution is used for time-independent variation}

\item{numApprox}{the number of normal distributions to approximate the Laplace distribution}

\item{margin}{the total probability mass that the number of jumps omitted in a compound Poisson process}

\item{numCores}{the number of cores to run in parallel}

\item{asymptotic}{the threshold of expected number of jumps on a branch beyond which normal distribution is assumed}
}
\value{
\verb{$summary} is a data frame listing the means and variances of the hidden states

\verb{$error} is a list of error distributions where each element is a data frame
}
\description{
Predict hidden states under the pulsed evolution model in a leave-one-out cross-validation
}
