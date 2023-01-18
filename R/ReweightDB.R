# @file ReweightDB.R
#
# Copyright 2021 KI Research Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @import CVXR
#' @import glue
NULL

#' @title Reweight an internal database to match the means of an external one.
#'
#' @description
#'
#' This function receives a data frame Z of an internal DB and a vector mu of
#' means from an external DBs. The elements of mu correspond to columns of Z.
#' It returns a set of sample weights such that the weighted means of the
#' columns of Z are as close as possible to the elements of mu while
#' minimizing the divergence between the distribution of the weights and the
#' uniform distribution.
#'
#' @param Z a data frame where every row stores a sample of the internal
#' databse.
#' @param mu a vector of means of the external dataset.
#' @param divergence 'entropy' or 'chi2'.
#' 'entropy' directs the algorithm to minimize the negative entropy,
#' \eqn{-\sum_i w_i \log w_i}.
#' 'chi2' is \eqn{\sum_i{w_i-\frac{1}{n}}**2}
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#' @param minW minimum weight
#' @param distance distance between means, either 'l1' or 'l2'
#' @param optimizationMethod primal or dual. Currently dual works only with entropy divergence
#' @param solver solver for convex problems. One of "ECOS" "ECOS_BB" "SCS" "OSQP"
#' @param verbose a Boolean flag for output messages
#'
#' @return
#' A vector of weights
#'
#' @export
reweightByMeans <- function(
  Z, mu, divergence = 'entropy', lambda=1e-6, minSd=1e-4, minW=1e-6, distance = 'l2', optimizationMethod = 'primal',
  solver = 'ECOS', verbose=FALSE)
{
  # TODO in the regularized case, hyper param may depend on the number of features/samples
  # Find optimal weights
  if (optimizationMethod == 'primal')
    w_hat <- primalReweightByMeans(Z, mu, divergence, lambda, minSd, minW, distance, solver, verbose)
  else
    w_hat <- dualReweightByMeans(Z, mu, lambda, minSd, minW, solver, verbose)

  min_w <- min(w_hat)
  mean_w <- mean(w_hat)
  w_summary <- sprintf('mean(w) = %.2f (sd = %.2f, min = %.2g, max = %.2g)', mean_w, sd(w_hat), min_w, max(w_hat))
  ParallelLogger::logInfo(w_summary)
  if ((!is.na(min_w)) & (min_w < 0)) {
    ParallelLogger::logWarn("Trimming negative weights to zero")
    n <- nrow(Z)
    w_hat[w_hat<0] <- minW/n
  }
  return (w_hat)
}



#' @title a reweighting object.
#'
#' @description
#'
#' @param divergence 'entropy' or 'chi2'.
#' 'entropy' directs the algorithm to minimize the negative entropy,
#' \eqn{-\sum_i w_i \log w_i}.
#' 'chi2' is \eqn{\sum_i{w_i-\frac{1}{n}}**2}
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#' @param optimizationMethod primal or dual. Currently dual works only with entropy divergence
#'
#' @return
#' An object of class cvxWeightOptimizer
#'
#' @export
cvxWeightOptimizer <- function(lambda=1e-1, minSd=1e-4, divergence = 'entropy', optimizationMethod = 'dual') {
  l <- list(
    shortName = 'CVX',
    lambda=lambda,
    minSd=minSd,
    divergence = divergence,
    optimizationMethod = optimizationMethod,
    optimize=optimizeWeightCVX
  )
  class(l) <- 'cvxWeightOptimizer'
  return(l)
}


optimizeWeightCVX <- function(wOptimizer, Z, mu) {
  w_hat <- reweightByMeans(
    Z, mu, divergence = wOptimizer$divergence, lambda = wOptimizer$lambda, minSd = wOptimizer$minSd, minW = 0,
    solver='ECOS', verbose = T, optimizationMethod = wOptimizer$optimizationMethod)
  r <- list(w_hat=w_hat)
  return(r)
}



#' Reweight an internal database to match the means of an external via solving a corresponding optimization problem.
#'
#' @description
#'
#' This function receives a data frame Z of an internal DB and a vector mu of
#' means from an external DBs. The elements of mu correspond to columns of Z.
#' It returns a set of sample weights such that the weighted means of the
#' columns of Z are as close as possible to the elements of mu while
#' minimizing the divergence between the distribution of the weights and the
#' uniform distribution.
#'
#'
#' @param Z a data frame where every row stores a sample of the internal
#' databse.
#' @param mu a vector of means of the external dataset.
#' @param divergence 'entropy' or 'chi2'.
#' 'entropy' directs the algorithm to minimize the negative entropy,
#' \eqn{-\sum_i w_i \log w_i}.
#' 'chi2' is \eqn{\sum_i{w_i-\frac{1}{n}}**2}
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#' @param minW mimimum weight
#' @param distance distance between means, either 'l1' or 'l2'
#' @param solver solver for convex problem
#' @param verbose a boolean flag for output messages
#'
#' @return
#' A vector of weights
#'
primalReweightByMeans <- function(Z, mu, divergence,lambda, minSd, minW, distance, solver, verbose) {
  normalized <- normalizeDataAndExpectations(Z, mu, minSd)
  n <- nrow(normalized$Z)
  w <- CVXR::Variable(n, 1)

  normalized$mu <- as.vector(normalized$mu)
  normalized$Z <- as.matrix(normalized$Z)

  if (divergence == 'entropy')    fDivergence <- -mean(entr(w))
  else if (divergence == 'chi2')  fDivergence <- CVXR::norm2(w-(1/n)) ** 2
  else {
    ParallelLogger::logError(glue("unsuported divergence type {divergence}"))
    return(rep(NaN, n))
  }

  if (distance == 'l2')       expectationsDistance <- CVXR::norm2(t(normalized$Z) %*% w - normalized$mu)
  else if (distance == 'l1')  expectationsDistance <- CVXR::norm1(t(normalized$Z) %*% w - normalized$mu)
  else {
    ParallelLogger::logError(glue("unsuported distance type {distance}"))
    return(rep(NaN, n))
  }

  if (lambda > 0) {
    if (verbose) cat(glue('Reweighting using {divergence}, {distance}, lambda = {lambda}, minW = {minW}'), '\n')
    objective <- CVXR::Minimize(expectationsDistance + lambda*fDivergence)
    constr <- list(w >= minW, sum(w) == 1)
  } else {
    if (verbose) cat(glue('Reweighting using {divergence}, hard expectation constraints, minW = {minW}'), '\n')
    objective <- CVXR::Minimize(fDivergence)
    constr <- list(w >= minW, sum(w) == 1, (t(normalized$Z) %*% w) == normalized$mu)
  }
  problem <- Problem(objective, constraints = constr)
  result <- solve(problem, solver = solver, verbose = FALSE)
  # The status of the solution can be "optimal", "optimal_inaccurate", "infeasible", "infeasible_inaccurate",
  # "unbounded", "unbounded_inaccurate", or "solver_error"
  if (result$status != 'optimal') {
    warning(glue('non-optimal results, returning NaNs, data size = {n} * {length(mu)}'))
    w_hat <- rep(NaN, n)
  } else {
    w_hat <- result$getValue(w) * n
  }
  return (w_hat)
}


#' Reweight an internal database to match the means of an external one using
#' a dual formulation of a weight optimization problem.
#'
#' @description
#'
#' This function recieves a data frame Z of an internal DB and a vector mu of
#' means from an external DBs. The elements of mu correspond to columns of Z.
#' It returns a set of sample weights such that the weighted means of the
#' columns of Z are as close as possible to the elements of mu while
#' minimizing the divergence between the distribution of the weights and the
#' uniform distribution.
#'
#' @param Z a data frame where every row stores a sample of the internal databse.
#' @param mu a vector of means of the internal dataset.
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#' @param minW mimimum weight
#' @param solver solver for convex problem
#' @param verbose a boolean flag for output messages
#'
#' @return
#' A vector of weights
#'
dualReweightByMeans <- function(Z, mu, lambda, minSd, minW, solver, verbose) {
  normalized <- normalizeDataAndExpectations(Z, mu, minSd)
  m <- ncol(normalized$Z)
  n <- nrow(normalized$Z)

  nu <- CVXR::Variable(m+1)
  C <- rbind(t(normalized$Z), matrix(1,1,n))
  d <- c(normalized$mu, 1)
  objective <- CVXR::Minimize(t(d) %*% nu + exp(-1) * sum_entries(exp(- t(C) %*% nu)))
  if (lambda==0)
    problem <- Problem(objective)
  else {
    constr <- list(CVXR::norm2(nu[1:m]) <= (1/lambda))
    problem <- Problem(objective, constraints = constr)
  }
  result <- solve(problem, solver=solver)

  if (result$status != 'optimal') {
    warning(glue('non-optimal results, returning NaNs, data size = {n} * {length(mu)}'))
    w_hat <- rep(NaN, n)
  } else {
      nu_hat <- result$getValue(nu)
      if (verbose)
        cat('nu:', nu_hat[1], ',' , nu_hat[2],  ',' , nu_hat[3], ', ... ,', nu_hat[m],  ',',nu_hat[m+1], '\n')
      if (!any(is.na(nu_hat))) {
          w_hat <- exp(-1- t(C) %*% nu_hat) * n
      } else {
        warning(glue('Optimization resulted in NaNs, data size = {n} * {length(mu)}'))
        w_hat <- rep(NaN, n)
      }
  }
  return(w_hat)
}


#' Normalize data and expectations
#'
#' @description
#'
#' Normalize data and expectations.
#'
#' @param Z data frame
#' @param mu vector of expectations
#' @param minSd float value of minimum standard deviation. Columns with a smaller sd will be removed
#'
#' @return
#' list
#'
#' @export
normalizeDataAndExpectations <- function(Z, mu, minSd) {
  # Preprocess and filter columns with low variance
  muZ <- colMeans(Z)
  sdZ <- apply(Z, 2, sd)
  useCols <- sdZ>minSd # remove columns with low variance
  if (sum(!useCols) > 0) {
    ParallelLogger::logInfo(
      sprintf("Removing %d columns with low variance (<%.2g). Please examing this columns.", sum(!useCols), minSd))
    Z <- Z[,names(Z)[useCols]]  # remove low variance columns todo: flag a positivity issue
    sdZ <- sdZ[useCols]
    muZ <- muZ[useCols]
    mu <- mu[useCols]
  }
  # Normalize
  p <- ncol(Z)
  for (i in 1:p) Z[,i] <- (Z[,i]-muZ[i])/sdZ[i]
  mu <- (mu-muZ)/sdZ
  return(list(Z=Z, mu=mu, muZ=muZ, sdZ=sdZ))
}


#' @title Compute Table 1 like transformation
#'
#' @description
#'
#' Compute features and squared numeric feature multiplied by outcome and 1-outcome.
#'
#' @param X data frame
#' @param outcomeBalance boolean specifying if to use outcome in balancing
#' @param outcomeCol string specifying the name of the outcome column
#'
#' @return
#' A data frame
#'
#' @export
computeTable1LikeTransformation <- function(X, outcomeBalance, outcomeCol='Y') {
  # Add squares of numeric features
  is_numeric <- sapply(X, function(x) length(unique(x))>2)  # TODO - consider a more elegant way
  for (s in names(is_numeric[is_numeric])) {
    sNew <- paste(s, '_Table1T_squared', sep='')
    X[[sNew]] <- (X[[s]]**2)
  }
  # Convert binary CVXR::Variables to numeric 0-1
  is_factor <- sapply(X, is.factor)
  X[is_factor] <- sapply(X[is_factor], as.numeric)
  X[is_factor] <- X[is_factor] - 1
  # Add interactions with outcome
  if (outcomeBalance) {
    Z = data.frame(row.names = row.names(X))
    Z[['Y']] <- X[[outcomeCol]]  # TODO decide if to maintain this
    x_names <- names(X)
    ext_x <- x_names[x_names != outcomeCol]
    ext_x_y1 <- paste(ext_x, '_Table1T_times_y1', sep="")
    for (i in 1:length(ext_x)) {
      Z[[ext_x_y1[i]]] <- X[[outcomeCol]] * X[[ext_x[i]]]
    }
    ext_x_y0 <- paste(ext_x, '_Table1T_times_y0', sep="")
    for (i in 1:length(ext_x)) {
      Z[[ext_x_y0[i]]] <- (1-X[[outcomeCol]]) * X[[ext_x[i]]]
    }
  } else {
    Z <- X
  }
  return (Z)
}


#' @title create reweight settings
#'
#' @param divergence distributional divergence
#' @param lambda regularization
#' @param minSd minimum standard deviation of features
#' @param minW minimum weight
#' @param distance distance between mean vectors
#' @param optimizationMethod primal or dual
#' @param solver solver
#'
#' @return
#' An object of class \code{reweightSettings}
#'
#' @export

createReweightSettings <- function(
    divergence = 'entropy',
    lambda=1e-6,
    minSd=1e-4,
    minW=1e-6,
    distance = 'l2',
    optimizationMethod = 'primal',
    solver = 'ECOS'
    )
{

  legalDivergences <- c('entropy', 'chi2')
  if (!(divergence %in% legalDivergences))
    stop(glue('unsupported divergence {divergence}. Should be in ({legalDivergences}).'))
  legalOptimizationMethods <- c('primal', 'dual')
  if (!(optimizationMethod %in% legalOptimizationMethods))
    stop(glue('unsupported optimization method {optimizationMethod}. Should be in ({legalOptimizationMethods}).'))
  reweightSettings <- list(
    divergence = divergence,
    lambda = lambda,
    minSd = minSd,
    minW = minW,
    distance = distance,
    optimizationMethod = optimizationMethod,
    solver = solver
  )
  class(reweightSettings) <- 'reweightSettings'
  return(reweightSettings)

}
