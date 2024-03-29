% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvxWeightOptimizer.R
\name{reweightByMeans}
\alias{reweightByMeans}
\title{Reweight an internal database to match the means of an external one.}
\usage{
reweightByMeans(
  Z,
  mu,
  divergence = "entropy",
  lambda = 1e-06,
  minSd = 1e-04,
  minW = 1e-06,
  distance = "l2",
  optimizationMethod = "primal",
  solver = "ECOS",
  verbose = FALSE
)
}
\arguments{
\item{Z}{a data frame where every row stores a sample of the internal
databse.}

\item{mu}{a vector of means of the external dataset.}

\item{divergence}{'entropy' or 'chi2'.
'entropy' directs the algorithm to minimize the negative entropy,
\eqn{-\sum_i w_i \log w_i}.
'chi2' is \eqn{\sum_i{w_i-\frac{1}{n}}**2}}

\item{lambda}{lambda - regularization parameter}

\item{minSd}{minimum variance for a columns to be included}

\item{minW}{minimum weight}

\item{distance}{distance between means, either 'l1' or 'l2'}

\item{optimizationMethod}{primal or dual. Currently dual works only with entropy divergence}

\item{solver}{solver for convex problems. One of "ECOS" "ECOS_BB" "SCS" "OSQP"}

\item{verbose}{a Boolean flag for output messages}
}
\value{
A vector of weights
}
\description{
This function receives a data frame Z of an internal DB and a vector mu of
means from an external DBs. The elements of mu correspond to columns of Z.
It returns a set of sample weights such that the weighted means of the
columns of Z are as close as possible to the elements of mu while
minimizing the divergence between the distribution of the weights and the
uniform distribution.
}
