#' Optimization using minimum squared error criterion with initial tuning of the
#' base learning rate.
#'
#' @param alphas
#' @param minSd
#' @param w0
#' @param nTuneIter
#' @param nIter
#' @param outputDir
#' @param improveTh
#' @param momentumC
#' @param approxUpdate
#' @param maxErr
#'
#' @return an object of class seTunedWeightOptimizer
#'
#' @export
seTunedWeightOptimizer <- function(
    alphas, minSd=1e-4, w0 = NULL,  nTuneIter=50, nIter=1000, outputDir=NULL, improveTh=1e-4, momentumC=0.9,
    approxUpdate=F, maxErr=0.01)
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
    maxErr = maxErr
  )
  class(l) <- 'seTunedWeightOptimizer'
  return(l)
}


setInitialW <- function(wOptimizer, results) {
  ParallelLogger::logInfo(glue('Setting initial w: sum={sum(results$w)}, sd={sd(results$w)}, n0- = {sum(results$w<=0)}'))
  wOptimizer$w0 <- results$w
  return(wOptimizer)
}



#' Optimization using minimum squered error criterion with initial tuning of the
#' base learning rate.
#'
#' @param wOptimizer
#' @param Z
#' @param mu
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
    maxErr = wOptimizer$maxErr
  )

  o <- rep(NA, length(wOptimizer$alphas))
  for (i in 1:length(o)) {
    cOptimizer <- tOptimizer
    cOptimizer$alpha <- wOptimizer$alphas[i]
    r <- cOptimizer$optimize(cOptimizer, Z, mu)
    gc()
    if (is.na(r$err))
      break
    o[i] <- r$err
  }
  runId <- gsub(':', '-', Sys.time())
  runPrefix = glue('optimizeSEWeightsTuned{runId}')
  # png(filename = file.path(wOptimizer$outputDir, glue('{runPrefix} epsilon.png')))
  # plot(wOptimizer$alphas, o, type='l', log='x', xlab = 'Alpha')
  # dev.off()


  k <- which.min(o)

  fOptimizer <- tOptimizer
  fOptimizer$nIter <- wOptimizer$nIter
  fOptimizer$alpha <- wOptimizer$alphas[k]
  r <- fOptimizer$optimize(fOptimizer, Z, mu)
  gc()
  return(r)
}



seOptimizer <- function(
    minSd=1e-4, w0 = NULL,  nIter=1000, outputDir=NULL, improveTh=1e-4, alpha=0.1, momentumC=0.9, approxUpdate=F,
    maxErr = 0.01)
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
    w0 = w0
  )
  class(l) <- 'seOptimizer'
  return(l)
}




optimizeSEWeights <- function(wOptimizer, Z, mu) {

  # approxUpdate = wOptimizer$approxUpdate

  alpha <- wOptimizer$alpha
  cat('alpha =',alpha, '\n')
  n <- nrow(Z)
  m <- ncol(Z)
  Y <- Z[['Y']]
  momentumC <- wOptimizer$momentumC
  v <- rep(0, n)

  normalized <- normalizeDataAndExpectations(Z, mu, wOptimizer$minSd)  # TODO
  cat('optimizeSEWeights\n')

  maxIter <- wOptimizer$nIter
  l <- matrix(nrow = maxIter, ncol = 3)  # Optimization log
  colnames(l) <- c('Primal objective', 'KL', 'Max abs Err')  #  'sum w',

  if (!is.null(wOptimizer$w0)) {
    w <- wOptimizer$w0/sum(wOptimizer$w0)
    ParallelLogger::logInfo(glue('Using initial w: sum={sum(w)}, sd={sd(w)}, n0- = {sum(w<=0)}'))
  }
  else
    w <- rep(1/n, n)  # Start with uniform weights
  muHat <- t(Z) %*% w  # estimated means
  rr <- muHat[ ,1] - mu
  err <- sum(rr**2)
  kl <- log(n) + t(w) %*% log(w)
  l[1, 1] <- err
  l[1, 2] <- kl
  l[1, 3] <- max(abs(rr))
  ParallelLogger::logInfo(glue('Starting with \terr {err}\tkl {kl}\tmax {l[1,3]}'))
  Z <- as.matrix(Z)
  g <- Z %*% rr  # * 2  #
  for (i in 2:maxIter) {

    g <- Z %*% rr  # * 2  #
    v <- momentumC * v + (1-momentumC) * g
    # if (approxUpdate) {
    #   w <- w * (1 - alpha * v)
    #   w[w <= 1e-9] <- 1e-9
    # }
    # else
    w <- w * exp(- alpha * v)
    w <- w/sum(w)

    muHat <- t(Z) %*% w
    rr <- muHat[, 1] - mu
    err <- sum(rr**2)
    maxAbsErr <- max(abs(rr))

    l[i, 1] <- err
    l[i, 3] <- maxAbsErr
    if (i%%100==0) {
      cat(i, err, kl, '\n')
      kl <- log(n) + t(w) %*% log(w)
      l[i,2] <- kl
    }
    if ((err < wOptimizer$improveTh) && (maxAbsErr < wOptimizer$maxErr))
      break
  }
  idx0 <- w<=1e-9  # TODO -
  n0 <- sum(idx0)
  if (n0>0) {
    ParallelLogger::logWarn(glue('{n0} weights = 0'))
    w[idx0] <- 1e-9
    w <- w/sum(w)
  }
  l[i, 2] <- log(n) + t(w) %*% log(w)
  l <- l[1:i, ]
  r <- list(w_hat=n*w, err=err, log=l)
  ParallelLogger::logInfo(glue('Finished at {i}\terr {err}\tkl {l[i,2]}\tmax {l[i,3]}'))
  return(r)
}


