#' Optimization using minimum squared error criterion with initial tuning of the
#' base learning rate.
#'
#' @param epsilon0 learning rates vector
#' @param nTuneIter number of iteration for finding the correct epsilon
#' @param minSd minimum standard deviation for removing features
#' @param nu0 initial value
#' @param nIter maximum number of iterations
#' @param outputDir output directory
#' @param improveTh threshold for stopping
#' @param momentumC momentum coefficient
#' @param batchSize batch size
#' @param normalizeMu force Nu to represent a distribution
#' @param polyBeta coefficient for learning rate decay
#' @param maxAbsNu maximum absolute value of nu to avoid numeric overflow
#' @param maxNuNorm  regulaziation parameter to limit max nu norm
#' @param nProbe number of iterations for showing intermediate resualts and assessments
#' @param averaging average the parameters over iterations
#'
#' @return an object of class seTunedWeightOptimizer
#'
#' @export
sgdTunedWeightOptimizer <- function(
    epsilon0s, nTuneIter=50, minSd=1e-4, nu0 = NULL,  nIter=1000, batchSize=20000,
    normalizeNu=T, momentumC = 0.9, outputDir=NULL, polyBeta=1, maxAbsNu=16, maxNuNorm=10, improveTh=1e-4, nProbe=50,
    averaging=F)
{
  l <- list(
    shortName = 'DUAL-SGD',
    epsilon0s = epsilon0s,
    nTuneIter = nTuneIter,
    minSd = minSd,
    nu0 = nu0,
    nIter = nIter,
    batchSize = batchSize,
    normalizeNu = normalizeNu,
    momentumC = momentumC,
    optimize = optimizeWeightSGDTuned,
    setInitialValue = setInitialNu,
    outputDir = outputDir,
    polyBeta = polyBeta,
    maxAbsNu = maxAbsNu,
    maxNuNorm = maxNuNorm,
    improveTh = improveTh,
    nProbe = nProbe,
    averaging = averaging
  )
  class(l) <- 'sgdTunedWeightOptimizer'
  return(l)
}


setInitialNu <- function(wOptimizer, results) {
  wOptimizer$nu0 <- results$nu
  return(wOptimizer)
}


optimizeWeightSGDTuned <- function(wOptimizer, Z, mu) {

  n <- nrow(Z)
  Y <- Z[['Y']]

  Z <- as.matrix(Z)  # TODO do it before

  ParallelLogger::logInfo(glue('n = {n}'))
  balanceTrain <- balanceClassesByUpsampling(Y)  # TODO

  tOptimizer <- getTunedOptimizer(wOptimizer, Z, mu, balanceTrain)

  maxTrials <- 5
  epsilon0Reduce <- 0.5
  # Try to optimize until success or maxTrials
  for (i in 1:maxTrials) {
    r <- tOptimizer$optimize(
      tOptimizer, Z, mu, balanceTrain, verbose=T,  #  cAll, extendedMu
      improveTh=wOptimizer$improveTh, nProbe=wOptimizer$nProbe, averaging=wOptimizer$averaging)
    gc()
    if (!any(is.na(r$nu)) && max(abs(r$nu))<=tOptimizer$maxAbsNu) { # TODO configure && r$dualValue<objectives[k]) Handle this better
      m <- ncol(Z)
      w_hat <- exp(-1 - Z %*% r$nu[1:m] - r$nu[m+1])

      if (!is.na(wOptimizer$batchSize) ) {
        idx1 <- Y==1
        w_hat[idx1] <- w_hat[idx1] * balanceTrain$nDup
      }
      idx0 <- w_hat == 0
      n0 <- sum(idx0)
      if (n0>0) {
        ParallelLogger::logWarn(glue('{n0} weights = 0'))
        w_hat[idx0] <- 1e-9
      }
      w_hat <- w_hat/sum(w_hat)
      w_hat <- w_hat*length(w_hat)
      r$w_hat <- w_hat
      return(r)

    }
    else {
      ParallelLogger::logInfo(glue('Reducing e0 from {tOptimizer$epsilon0} by {epsilon0Reduce}'))
      tOptimizer$epsilon0 <- tOptimizer$epsilon0 * epsilon0Reduce
    }
  }
  return(r)
}


#' Generate information for class balancing during training and weight coputations
#'
#' @param y vector of outcomes
#'
#' @return a named list with \code{sampleIdx} and \code{nDup}
#'
balanceClassesByUpsampling <- function(y) {
  idx1 <- which(y==1)
  idx0 <- which(y==0)
  n0 <- length(idx0)
  n1 <- length(idx1)
  ParallelLogger::logInfo(glue('n1={n1}, n0={n0}'))
  samplesIdx <- 1:length(y)
  if (n1<n0/2) {  # TODO add a flag here
    nDup <- floor(n0/n1)
    ParallelLogger::logInfo(glue('duplicating class 1 by {nDup}'))
    for (i in 1:(nDup-1))
      samplesIdx <- c(samplesIdx, idx1)
  } else {
    nDup <- 1
  }
  n <- length(samplesIdx)
  shuffleIdx <- sample(n, n, replace = F)
  samplesIdx <- samplesIdx[shuffleIdx]
  ParallelLogger::logInfo(glue('Balance-augmented set size = {n}'))
  return(list(samplesIdx=samplesIdx, nDup=nDup))
}


exponentialEpsilonSchedule <- function(wOptimizer, n) {
  i <- 1:wOptimizer$nIter
  gamma <- 0.001
  epsilon <- wOptimizer$epsilon0 * exp(-gamma*i)
  return(epsilon)
}


polynomialEpsilonSchedule <- function(wOptimizer, n) {
  i <- 1:wOptimizer$nIter
  alpha = 0.5
  beta = wOptimizer$polyBeta
  epsilon <- wOptimizer$epsilon0*(beta * i + 1)**(-alpha)
}

#' Get a tuned optimizer
#'
#' @description select a step size that gives the best objective after a few
#' starting iterations
#'
#' @param wOptimizer optimization object of type \code{sgdTunedWeightOptimizer}
#' @param Z feature matrix
#' @param mu statistics vector
#' @param balanceTrain balancing object
#'
#' @return tuned optimizer
#'
getTunedOptimizer <- function(wOptimizer, Z, mu, balanceTrain) {

  runId <- gsub(':', '-', Sys.time())
  runPrefix = glue('optimizeWeightSGDTuned{runId}')
  m <- length(mu)

  swOptimizer <- sgdWeightOptimizer(
    minSd = wOptimizer$minSd,
    nu0 = wOptimizer$nu0,
    nIter = NaN,  #
    epsilon0 = NaN, #
    polyBeta = wOptimizer$polyBeta,
    batchSize = wOptimizer$batchSize,
    normalizeNu = wOptimizer$normalizeNu,
    momentumC = wOptimizer$momentumC,
    maxAbsNu = wOptimizer$maxAbsNu,
    maxNuNorm = wOptimizer$maxNuNorm
  )

  nTune <- length(wOptimizer$epsilon0s)
  objectives <- rep(NaN, nTune)
  epsilon0s <- sort(wOptimizer$epsilon0s)
  for (i in 1:nTune) {
    iOptimizer <- swOptimizer
    iOptimizer$nIter <- wOptimizer$nTuneIter  #
    iOptimizer$epsilon0 <- epsilon0s[i] #
    r  <- iOptimizer$optimize(
      iOptimizer, Z, mu, balanceTrain, improveTh=wOptimizer$improveTh, averaging=wOptimizer$averaging)
    gc()
    if (any(is.na(r$nu)) || max(abs(r$nu))>wOptimizer$maxAbsNu)  # TODO assuming that if epsilon0 overflows so will larger values
      break
    r$nu[m+1, 1] <- -1 + log(sum(exp(-t(Z[1:m,]) %*% r$nu[1:m,1])))  # Normalize
    objectives[i] <- r$log[nrow(r$log), 1]  # TODO verify that this is not batch dependant
  }
  # png(filename = file.path(wOptimizer$outputDir, glue('{runPrefix} epsilon.png')))
  # plot(epsilon0s, objectives, type='l', log='x')
  # dev.off()
  k <- which.min(objectives)
  #if (k>1)
  #   k <- k-1  # move one to the left to make the choice more stable
  ParallelLogger::logInfo(glue('epsilon0 = {epsilon0s[k]}'))
  tOptimizer <- swOptimizer
  tOptimizer$nIter = wOptimizer$nIter  #
  tOptimizer$epsilon0 = epsilon0s[k]  #
  return(tOptimizer)

}


#' Optimization using stochastic gradient descent
#'
#' @param wOptimizer and optimizer object
#' @param Z feature matrix
#' @param muExt external statistics
#' @param balanceTrain balancing object
#' @param verbose
#'
#' @return optimization results. A list
#'
#'
optimizeWeightSGD <- function(
    wOptimizer, Z, muExt, balanceTrain, improveTh = 1e-4, nProbe = 50, averaging=F) {

  # Init params
  n <- nrow(Z)
  m <- ncol(Z)
  batchSize <- min(wOptimizer$batchSize, n)
  nUpsampled <- length(balanceTrain$samplesIdx)
  nEpoches <- (batchSize * wOptimizer$nIter) / nUpsampled
  batchesPerEpoch <- floor(nUpsampled/batchSize)
  epsilon <- polynomialEpsilonSchedule(wOptimizer, n)   # or exponentialEpsilonSchedule
  l <- matrix(nrow = wOptimizer$nIter, ncol=4)  # Optimization log
  colnames(l) <- c('Dual objective', 'Primal objective', 'Norm', 'Other')  #  'sum w',

  # Init nu
  nu <- wOptimizer$nu0
  if (is.null(nu)) {
    # Initialize to represent uniform weights
    nu <- matrix(data=rep(0, m+1))
    nu[m+1, 1] <- log(batchSize) - 1 # Normalize ... or n
  }

  # Init momentum. TODO consider optimizing the gradient
  v <- rep(0, m+1)
  nuBar <- rep(0, m+1)  # TODO this is for iterative averaging
  nRescales <- 0
  g <- rep(NA, m+1)

  for (i in 1:wOptimizer$nIter) {

    # Epoch-wise monitoring and minibatch selection
    if (!is.na(batchSize)) {  # TODO check if this should be n
      batchIdx <- (i-1)%%batchesPerEpoch
      if (batchIdx==0) {
        shuffled <- balanceTrain$samplesIdx[sample(nUpsampled, nUpsampled, replace = F)]  # TODO
      }
      idx <- shuffled[(batchSize*batchIdx+1):(batchSize*(batchIdx+1))]
      Zi <- Z[idx, ]  # Ci <- C[, idx]
    } else
      Zi <- Z #  Ci <- C

    # Step
    nu <- nu + v  # option 1
    normNu <- norm(nu[1:m], type='2')  # Check norm and rescale if indicated
    if (normNu > wOptimizer$maxNuNorm) {
      nRescales = nRescales + 1
      nu[1:m] <- nu[1:m]*wOptimizer$maxNuNorm/normNu
    }

    nuBar <- 1/i * (nu) + (1-1/i) * nuBar  # TODO for 'averaging' consider this

    # wi <- exp(-t(Ci[1:m, ]) %*% nu[1:m,1])
    wi <- exp(-Zi %*% nu[1:m,1])
    if (wOptimizer$normalizeNu)
      nu[m+1, 1] <- -1 + log(sum(wi))
    if (any(is.na(nu)) || max(abs(nu))>wOptimizer$maxAbsNu) {
      ParallelLogger::logInfo(glue('Optimization overflow nu = {nu[m+1]}'))
      return(list(log=l, nu=nu))
    }
    wi <- wi*exp(-1-nu[m+1, 1])
    res <- (t(Zi) %*% wi)-muExt  # ...
    # g <- extendedMu - rowSums(Ci %*% wi)  # Dual gradient same as dualGradient(nu, Ci, extendedMu)
    g[1:m] <- muExt - rowSums(t(Zi) %*% wi)
    g[m+1] <- 1 - sum(wi)

    v <- wOptimizer$momentumC * v - epsilon[i] * g
    l[i, 1] <- t(muExt) %*% nu[1:m] + nu[m+1] + sum(wi)  # dual objective
    l[i, 2] <- sum(res**2)   # o['risk']  # TODO Notice this should we normalize?
    l[i, 3] <- normNu
    l[i, 4] <- max(abs(res))

    if (l[i, 'Primal objective'] < improveTh)
      break
  }
  l <- l[1:i, ]

  if (nRescales>0)
    ParallelLogger::logInfo(glue('Rescaled nu to {wOptimizer$maxNuNorm} {nRescales} times'))
  ParallelLogger::logInfo(glue('Finished at {i}\tdual {l[i,1]}\tprimal {l[i,2]}\tnorm {l[i,3]}\tmax {l[i,4]}'))

  if (averaging)
    return(list(log=l, nu=nuBar, err=l[i,2]))
  else
    return(list(log=l, nu=nu, err=l[i,2]))
}


#' 'Base class' of sgdWeightOptimizer
#'
stochasticWeightOptimizer <- function(minSd=1e-4, nu0 = NaN,  nIter=1000, epsilon0=1e-1, batchSize=1e4) {
  l <- list(
    minSd = minSd,
    nu0 = nu0,
    nIter = nIter,
    epsilon0 = epsilon0,
    batchSize = batchSize
  )
  return(l)
}


sgdWeightOptimizer <- function(
    minSd=1e-4, nu0 = NaN,  nIter=1000, epsilon0=1e-1, batchSize=20000,
    normalizeNu=F, momentumC = 0.9, outputDir = NULL, polyBeta=1, maxAbsNu=20, maxNuNorm=10)
{
  l <- stochasticWeightOptimizer(
    minSd = minSd,
    nu0 = nu0,
    nIter = nIter,
    epsilon0 = epsilon0,
    batchSize = batchSize
  )
  l$momentumC <- momentumC
  l$normalizeNu <- normalizeNu
  l$optimize <- optimizeWeightSGD
  l$outputDir <- outputDir
  l$polyBeta <- polyBeta
  l$maxAbsNu <- maxAbsNu
  l$maxNuNorm <- maxNuNorm

  class(l) <- 'sgdWeightOptimizer'
  return(l)
}
