% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RasperGade_diagnosis.R
\name{analyzeResidualErrorByPPplot}
\alias{analyzeResidualErrorByPPplot}
\alias{analyzeResidualErrorByHeteroscedasticity}
\alias{analyzeResidualErrorByCI}
\title{Analyze quality of uncertainty in ancestral or hidden state prediction}
\usage{
analyzeResidualErrorByPPplot(
  trait,
  pred,
  error,
  epsilon = 0,
  laplace = FALSE,
  discrete = FALSE,
  n = 1000
)

analyzeResidualErrorByHeteroscedasticity(
  trait,
  pred,
  error,
  epsilon = 0,
  laplace = FALSE,
  bin = 100,
  n = 1000
)

analyzeResidualErrorByCI(
  trait,
  pred,
  error,
  epsilon = 0,
  laplace = FALSE,
  bin = 20,
  alpha = 0.05
)
}
\arguments{
\item{trait}{true trait values}

\item{pred}{point estimates of ancestral or hidden states}

\item{error}{uncertainty of predictions}

\item{epsilon}{variance of time-independent variation}

\item{laplace}{logical, if true, time-independent variation follows a Laplace distribution}

\item{discrete}{logical, if true, results should be analyzed as integers}

\item{bin}{number of bins in the test}
}
\description{
Analyze quality of uncertainty in ancestral or hidden state prediction from different perspectives
}
