% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preDiagnostics.R
\name{getHighlySkewedBinaryFeatures}
\alias{getHighlySkewedBinaryFeatures}
\title{Get highly skewed binary features}
\usage{
getHighlySkewedBinaryFeatures(
  mu,
  z,
  minNumReport = 20,
  maxDiff = 0.01,
  maxSubset = 20000
)
}
\arguments{
\item{mu}{expectations}

\item{z}{data}

\item{minNumReport}{minimum number of samples to report in the logger}

\item{maxDiff}{max difference in the unary case}

\item{maxSubset}{assumption about the maximum subset size

TODO better treatment cases in which binary values may not be 0 or 1. Transform to this form if needed}
}
\value{
a list
}
\description{
Get highly skewed binary features
}
