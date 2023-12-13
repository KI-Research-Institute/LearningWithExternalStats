#' An optimizer that uses minimum squared error criterion with initial tuning of the
#' base learning rate.
#'
#' @param outcomeCol name of outcome column
#' @param alphas a set of candidate baseline rates
#' @param minSd minimum standard deviation of feature
#' @param w0 initial weights vector, if NULL the weights will be initialized to uniform
#' @param nTuneIter number of iterations for tuning
#' @param nIter maximum number of iterations
#' @param outputDir output directory for logging
#' @param absTol required accuracy
#' @param momentumC momentum parameter
#' @param approxUpdate using 1+x instead of exp(x), currently disabled
#' @param absMaxUnivariateTol alternative stopping criterion: error of the maximum statistic
#' @param experimental use experimental method (Boolean)
#' @param maxSuccessMSE maximum MSE that is considered successful convergence
#'
#' @return an object of class \code{seTunedWeightOptimizer}
#'
#' @export
seTunedWeightOptimizer <- function(
    outcomeCol,
    alphas = c(0.01, 0.03, 0.1, 0.3, 0.5), minSd=1e-4, w0 = NULL,  nTuneIter=50, nIter=2000, outputDir=NULL,
    absTol=1e-8, momentumC=0.9, approxUpdate=F, absMaxUnivariateTol=1e-9, experimental=F, maxSuccessMSE=1e-5)
{
  l <- list(
    outcomeCol = outcomeCol,
    shortName = 'W-MSE',
    alphas = alphas,
    minSd = minSd,
    w0 = w0,
    nTuneIter = nTuneIter,
    nIter = nIter,
    outputDir = outputDir,
    absTol = absTol,
    alpha = NULL,
    momentumC = momentumC,
    approxUpdate = approxUpdate,
    optimize = optimizeSEWeightsTuned,
    setInitialValue = setInitialW,
    absMaxUnivariateTol = absMaxUnivariateTol,
    experimental = experimental,
    maxSuccessMSE = maxSuccessMSE
  )
  class(l) <- 'seTunedWeightOptimizer'
  return(l)
}


setInitialW <- function(wOptimizer, results) {
  # ParallelLogger::logInfo(glue('Setting initial w: sum={sum(results$w)}, sd={sd(results$w)}, n0- = {sum(results$w<=0)}'))
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
  # Normalize if not all columns are binary
  isBinary <- apply(Z, 2,function(x) {all(x %in% c(0, 1))})
  if (!all(isBinary)) {
    normalized <- normalizeDataAndExpectations(Z, mu, wOptimizer$minSd)  # TODO
    rm(list=c('Z', 'mu'))
    gc()
    Z <- normalized$Z
    mu <- normalized$mu
  }
  Z <- as.matrix(Z)

  # Get initinal error
  rr <- colMeans(Z) - mu
  err <- sum(rr**2)
  # ParallelLogger::logInfo(glue('Initial err {err}'))

  # Initialize optimizer
  tOptimizer <- seOptimizer(
    minSd = wOptimizer$minSd,
    nIter = wOptimizer$nTuneIter,
    outputDir = wOptimizer$outputDir,
    absTol = wOptimizer$absTol,
    w0 = wOptimizer$w0,
    alpha = NULL,
    momentumC = wOptimizer$momentumC,
    approxUpdate = wOptimizer$approxUpdate,
    absMaxUnivariateTol = wOptimizer$absMaxUnivariateTol,
    experimental = wOptimizer$experimental,
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
    if (is.null(r$err) || is.null(r$maxAbsErr) || is.na(r$err) || is.na(r$maxAbsErr)) {
      ParallelLogger::logInfo(glue('No results for alpha={cOptimizer$alpha}, skipping'))
      break
    }
    if ((r$err < wOptimizer$absTol) || (r$maxAbsErr < wOptimizer$absMaxUnivariateTol)) {
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
  fOptimizer$absTol <- wOptimizer$absTol
  fOptimizer$nIter <- wOptimizer$nIter
  fOptimizer$alpha <- wOptimizer$alphas[k]

  w_hat <- tuneResults[[k]]$w_hat
  fOptimizer$w0 <- w_hat/sum(w_hat)

  r <- fOptimizer$optimize(fOptimizer, Z, mu)
  gc()
  l <- r$log
  i <- nrow(l)
  totalIter <- totalIter + i
  r$totalIter <- totalIter
  ParallelLogger::logInfo(glue('Finished at {totalIter}\terr {r$err}\tkl {l[i,2]}\tmax {l[i,3]} status={r$status}'))
  return(r)
}



seOptimizer <- function(
    minSd=1e-4, w0 = NULL,  nIter=1000, outputDir=NULL, absTol=1e-9, alpha=0.1, momentumC=0.9, approxUpdate=F,
    absMaxUnivariateTol = 1e-5, experimental=F, maxSuccessMSE=0.1)
{
  l <- list(
    minSd = minSd,
    nIter = nIter,
    optimize = optimizeSEWeights,
    outputDir = outputDir,
    absTol = absTol,
    alpha = alpha,
    momentumC=momentumC,
    approxUpdate = approxUpdate,
    absMaxUnivariateTol = absMaxUnivariateTol,
    w0 = w0,
    experimental = experimental,
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
  Y <- Z[ ,wOptimizer$outcomeCol]
  momentumC <- wOptimizer$momentumC
  v <- rep(0, n)  # momentum velocity

  maxIter <- wOptimizer$nIter
  l <- matrix(nrow = maxIter, ncol = 3)  # Optimization log
  colnames(l) <- c('Primal objective', 'KL', 'Max abs Err')  #  'sum w',

  if (!is.null(wOptimizer$w0)) {
    w <- wOptimizer$w0/sum(wOptimizer$w0)
  }
  else
    w <- rep(1/n, n)  # Start with uniform weights

  # a heuristic initialization that sets weights according to Y
  idx1 <- Y==1
  idx0 <- Y==0
  n1 <- sum(idx1)
  n0 <- sum(idx0)
  w[idx1] <- (w[idx1]/sum(w[idx1])) * mu[wOptimizer$outcomeCol]  # TODO this relies on an assmptin the mu was not rescaled
  w[idx0] <- (w[idx0]/sum(w[idx0])) * (1-mu[wOptimizer$outcomeCol])
  # ParallelLogger::logInfo(glue('Using initial w: sum={sum(w)}, sd={sd(w)}, n0- = {sum(w<=0)}'))
  # First iteration information
  muHat <- t(Z) %*% w  # estimated means
  rr <- muHat[ ,1] - mu
  err <- sum(rr**2)

  kl <- log(n) + t(w) %*% log(w)  # TODO change the reference to a class wise
  l[1, 1] <- err
  l[1, 2] <- kl
  l[1, 3] <- max(abs(rr))
  # ParallelLogger::logInfo(glue('Starting with \terr {err}\tkl {kl}\tmax {l[1,3]}'))
  g <- Z %*% rr  # * 2  #
  for (i in 2:maxIter) {

    g <- Z %*% rr  # ignoring factor 2
    v <- momentumC * v + g
    w <- w * exp(- alpha * v)
    w <- w/sum(w)

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
      ParallelLogger::logWarn('NA objective')
      return(list(w_hat=w*n, err=Inf, log=l, status = 'NA-objective'))
    }
    if ((err < wOptimizer$absTol) || (maxAbsErr < wOptimizer$absMaxUnivariateTol))
      break
  }
  minWTh <- 1e-10
  idx0 <- w<=minWTh  # TODO -
  n0 <- sum(idx0)
  if (n0>0) {
    ParallelLogger::logInfo(glue('{n0} weights = 0'))
    w[idx0] <- minWTh
    w <- w/sum(w)
  }
  l[i, 2] <- log(n) + t(w) %*% log(w)
  l <- l[1:i, ]

  if (err<wOptimizer$maxSuccessMSE)
    status = 'Success'
  else
    status = 'Not-converged'
  r <- list(w_hat=n*w, err=err, log=l, status=status, maxAbsErr=maxAbsErr)
  return(r)
}


