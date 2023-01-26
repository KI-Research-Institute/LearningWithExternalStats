#' An optimizer that uses minimum squared error criterion with initial tuning of the
#' base learning rate.
#'
#' @param alphas a set of candidate baseline rates
#' @param minSd minimum standard deviation of feature
#' @param w0 initial weights vector, if NULL the weights will be initialized to uniform
#' @param nTuneIter number of iterations for tuning
#' @param nIter maximum number of iterations
#' @param outputDir output directory for logging
#' @param improveTh required accuracy
#' @param momentumC momentum parameter
#' @param approxUpdate using 1+x instead of exp(x), currently disabled
#' @param maxErr alternative stopping criterion: error of the maximum statistic
#' @param nesterov use Nesterov momentum (Boolean)
#' @param maxSuccessMSE maximum MSE that is considered successful convergence
#'
#' @return an object of class \code{seTunedWeightOptimizer}
#'
#' @export
seTunedWeightOptimizer <- function(
    alphas = c(0.01, 0.03, 0.1, 0.3, 0.5), minSd=1e-4, w0 = NULL,  nTuneIter=25, nIter=500, outputDir=NULL,
    improveTh=1e-4, momentumC=0.9, approxUpdate=F, maxErr=1, nesterov=F, maxSuccessMSE=0.1)
{
  l <- list(
    shortName = 'W-MSE',
    alphas = alphas,
    minSd = minSd,
    w0 = w0,
    nTuneIter = nTuneIter,
    nIter = nIter,
    outputDir = outputDir,
    improveTh = improveTh,
    alpha = NULL,
    momentumC = momentumC,
    approxUpdate = approxUpdate,
    optimize = optimizeSEWeightsTuned,
    setInitialValue = setInitialW,
    maxErr = maxErr,
    nesterov = nesterov,
    maxSuccessMSE = maxSuccessMSE
  )
  class(l) <- 'seTunedWeightOptimizer'
  return(l)
}


setInitialW <- function(wOptimizer, results) {
  ParallelLogger::logInfo(glue('Setting initial w: sum={sum(results$w)}, sd={sd(results$w)}, n0- = {sum(results$w<=0)}'))
  wOptimizer$w0 <- results$w
  return(wOptimizer)
}



#' Optimization using minimum squeread error criterion with initial tuning of the
#' base learning rate.
#'
#' @param wOptimizer optimizer object of type \code{seTunedWeightOptimizer}
#' @param Z feature matrix
#' @param mu external expectations vector
#'
#' @return optimization results. A list
#'
optimizeSEWeightsTuned <- function(wOptimizer, Z, mu)
{
  tOptimizer <- seOptimizer(
    minSd = wOptimizer$minSd,
    nIter = wOptimizer$nTuneIter,
    outputDir = wOptimizer$outputDir,
    improveTh = wOptimizer$improveTh,
    w0 = wOptimizer$w0,
    alpha = NULL,
    momentumC = wOptimizer$momentumC,
    approxUpdate = wOptimizer$approxUpdate,
    maxErr = wOptimizer$maxErr,
    nesterov = wOptimizer$nesterov,
    maxSuccessMSE = wOptimizer$maxSuccessMSE
  )

  o <- rep(NA, length(wOptimizer$alphas))
  tuneResults <- vector(mode = 'list', length = length(wOptimizer$alphas))
  totalIter <- 0
  for (i in 1:length(o)) {
    cOptimizer <- tOptimizer
    cOptimizer$alpha <- wOptimizer$alphas[i]
    r <- cOptimizer$optimize(cOptimizer, Z, mu)
    gc()
    totalIter <- totalIter + nrow(r$log)
    if (is.na(r$err))
      break
    if (r$err < wOptimizer$improveTh) {
      r$totalIter <- totalIter
      return(r)
    }
    o[i] <- r$err
    tuneResults[[i]] <- r
  }

  if (F) {  # For tuning
    runId <- gsub(':', '-', Sys.time())
    runPrefix = glue('optimizeSEWeightsTuned{runId}')
    png(filename = file.path(wOptimizer$outputDir, glue('{runPrefix} epsilon.png')))
    plot(wOptimizer$alphas, o, type='l', log='x', xlab = 'Alpha')
    dev.off()
  }

  k <- which.min(o)

  fOptimizer <- tOptimizer
  fOptimizer$nIter <- wOptimizer$nIter
  fOptimizer$alpha <- wOptimizer$alphas[k]

  w_hat <- tuneResults[[k]]$w_hat
  fOptimizer$w0 <- w_hat/sum(w_hat)

  r <- fOptimizer$optimize(fOptimizer, Z, mu)
  gc()
  totalIter <- totalIter + nrow(r$log)
  r$totalIter <- totalIter
  return(r)
}



seOptimizer <- function(
    minSd=1e-4, w0 = NULL,  nIter=1000, outputDir=NULL, improveTh=1e-4, alpha=0.1, momentumC=0.9, approxUpdate=F,
    maxErr = 0.01, nesterov=F, maxSuccessMSE=0.1)
{
  l <- list(
    minSd = minSd,
    nIter = nIter,
    optimize = optimizeSEWeights,
    outputDir = outputDir,
    improveTh = improveTh,
    alpha = alpha,
    momentumC=momentumC,
    approxUpdate = approxUpdate,
    maxErr = maxErr,
    w0 = w0,
    nesterov = nesterov,
    maxSuccessMSE = maxSuccessMSE
  )
  class(l) <- 'seOptimizer'
  return(l)
}




optimizeSEWeights <- function(wOptimizer, Z, mu) {

  # TODO deprecate approxUpdate = wOptimizer$approxUpdate

  alpha <- wOptimizer$alpha
  n <- nrow(Z)
  m <- ncol(Z)
  Y <- Z[['Y']]
  momentumC <- wOptimizer$momentumC
  v <- rep(0, n)  # momentum velocity

  # normalized <- normalizeDataAndExpectations(Z, mu, wOptimizer$minSd)  # TODO

  maxIter <- wOptimizer$nIter
  l <- matrix(nrow = maxIter, ncol = 3)  # Optimization log
  colnames(l) <- c('Primal objective', 'KL', 'Max abs Err')  #  'sum w',

  if (!is.null(wOptimizer$w0)) {
    w <- wOptimizer$w0/sum(wOptimizer$w0)
  }
  else
    w <- rep(1/n, n)  # Start with uniform weights

  Z <- as.matrix(Z)
  # a heuristic initialization that sets weights according to Y
  idx1 <- Y==1
  idx0 <- Y==0
  n1 <- sum(idx1)
  n0 <- sum(idx0)
  w[idx1] <- (w[idx1]/sum(w[idx1])) * mu['Y']  # TODO this relies on an assmptin the mu was not rescaled
  w[idx0] <- (w[idx0]/sum(w[idx0])) * (1-mu['Y'])
  # ParallelLogger::logInfo(glue('Using initial w: sum={sum(w)}, sd={sd(w)}, n0- = {sum(w<=0)}'))
  # First iteration information
  muHat <- t(Z) %*% w  # estimated means
  rr <- muHat[ ,1] - mu
  err <- sum(rr**2)
  kl <- log(n) + t(w) %*% log(w)  # TODO change the reference to a class wise
  l[1, 1] <- err
  l[1, 2] <- kl
  l[1, 3] <- max(abs(rr))
  ParallelLogger::logInfo(glue('Starting with \terr {err}\tkl {kl}\tmax {l[1,3]}'))
  g <- Z %*% rr  # * 2  #
  for (i in 2:maxIter) {

    if (!wOptimizer$nesterov) {
      g <- Z %*% rr  # ignoring factor 2
      v <- momentumC * v + g
      w <- w * exp(- alpha * v)
      w <- w/sum(w)

    } else {  # Nesterov momentum
      # First advance by the momentum and then compute the gradient there
      wbar <- w * exp(momentumC * v)
      wbar <- wbar/sum(wbar)
      muHat <- t(Z) %*% wbar
      rr <- muHat[, 1] - mu
      g <- Z %*% rr  # * 2  #

      v <- momentumC * v - alpha * g
      w <- w * exp(v)
      w <- w/sum(w)
    }
    muHat <- t(Z) %*% w
    rr <- muHat[, 1] - mu
    err <- sum(rr**2)
    maxAbsErr <- max(abs(rr))

    l[i, 1] <- err
    l[i, 3] <- maxAbsErr
    if (i%%100==0) {
      kl <- log(n) + t(w) %*% log(w)
      l[i,2] <- kl
    }
    if (is.na(err)) {
      return(list(w_hat=w*n, err=Inf, log=l))
    }
    if ((err < wOptimizer$improveTh) && (maxAbsErr < wOptimizer$maxErr))
      break
  }
  idx0 <- w<=1e-9  # TODO -
  n0 <- sum(idx0)
  if (n0>0) {
    ParallelLogger::logInfo(glue('{n0} weights = 0'))
    w[idx0] <- 1e-9
    w <- w/sum(w)
  }
  l[i, 2] <- log(n) + t(w) %*% log(w)
  l <- l[1:i, ]

  if (err<wOptimizer$maxSuccessMSE)
    status = 'Success'
  else
    status = 'Not-converged'
  r <- list(w_hat=n*w, err=err, log=l, status=status)
  ParallelLogger::logInfo(glue('Finished at {i}\terr {err}\tkl {l[i,2]}\tmax {l[i,3]} status={status}'))
  return(r)
}


