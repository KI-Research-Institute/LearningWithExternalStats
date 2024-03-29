% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mseWeightOptimizer.R
\name{seTunedWeightOptimizer}
\alias{seTunedWeightOptimizer}
\title{An optimizer that uses minimum squared error criterion with initial tuning of the
base learning rate.}
\usage{
seTunedWeightOptimizer(
  outcomeCol,
  alphas = c(0.01, 0.03, 0.1, 0.3, 0.5),
  minSd = 1e-04,
  w0 = NULL,
  nTuneIter = 50,
  nIter = 2000,
  outputDir = NULL,
  absTol = 1e-08,
  momentumC = 0.9,
  approxUpdate = F,
  absMaxUnivariateTol = 1e-09,
  experimental = F,
  maxSuccessMSE = 1e-05
)
}
\arguments{
\item{outcomeCol}{name of outcome column}

\item{alphas}{a set of candidate baseline rates}

\item{minSd}{minimum standard deviation of feature}

\item{w0}{initial weights vector, if NULL the weights will be initialized to uniform}

\item{nTuneIter}{number of iterations for tuning}

\item{nIter}{maximum number of iterations}

\item{outputDir}{output directory for logging}

\item{absTol}{required accuracy}

\item{momentumC}{momentum parameter}

\item{approxUpdate}{using 1+x instead of exp(x), currently disabled}

\item{absMaxUnivariateTol}{alternative stopping criterion: error of the maximum statistic}

\item{experimental}{use experimental method (Boolean)}

\item{maxSuccessMSE}{maximum MSE that is considered successful convergence}
}
\value{
an object of class \code{seTunedWeightOptimizer}
}
\description{
An optimizer that uses minimum squared error criterion with initial tuning of the
base learning rate.
}
