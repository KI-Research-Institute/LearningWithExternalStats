#' @importFrom stats sd quantile
#' @importFrom pROC auc roc
#' @importFrom WeightedROC WeightedAUC WeightedROC
#' @import ParallelLogger
NULL


#' evaluation from external stats

#' @title Create settings for \code{estimateExternalPerformanceFromStatistics}
#'
#' @param divergence 'entropy' or 'chi2'.
#' 'entropy' directs the algorithm to minimize the negative entropy, \eqn{\sum_i w_i \log w_i}.
#' 'chi2' is \eqn{\sum_i(w_i-\frac{1}{n})^2}
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#' @param minW minimum weight
#' @param distance distance between mean vectors, l1 or l2
#' @param optimizationMethod primal or dual. Currently dual works only with entropy divergence
#' @param solver solver ("ECOS" "ECOS_BB" "SCS" "OSQP")
#' @param nMaxReweight maximum number of samples for re-weighting. If this number is smaller than the number of
#' observations, then the estimator acts on sub-samples of the data-set.
#' @param nRepetitions number of repetitions with sub-samples.
#' @param stratifiedSampling attempt to sample a balanced number of samples for both outcome groups.
#' @param maxProp maximum proportion between external and internal means.
#' @param maxDiff maximum difference between  external and internal fequencies in case of unary variables or
#' binary variables with high proportion.
#' @param maxWSMD maximum allowed weighted standardized mean difference
#' @param outputDir output directory for logging
#' @param maxCores maximum number of cores for parallel processing of multiple sub-samples.
#'
#' @return
#' An object of class \code{externalEstimatorSettings}
#'
#' @export
createExternalEstimatorSettings <- function(
    divergence = 'entropy',
    lambda=1e-6,
    minSd=1e-4,
    minW=1e-6,
    distance = 'l2',
    optimizationMethod = 'primal',
    solver = 'ECOS',
    nMaxReweight = 10000,
    nRepetitions = 1,
    stratifiedSampling = T,
    maxProp = 500,
    maxDiff = 0.01,
    maxWSMD = 0.2,
    outputDir = getwd(),
    maxCores = 1
)
{
  externalEstimatorSettings <- list(
    divergence = divergence,
    lambda = lambda,
    minSd = minSd,
    minW = minW,
    distance = distance,
    optimizationMethod = optimizationMethod,
    solver = solver,
    nMaxReweight = nMaxReweight,
    nRepetitions = nRepetitions,
    stratifiedSampling = stratifiedSampling,
    maxProp = maxProp,
    maxDiff = maxDiff,
    maxWSMD = maxWSMD,
    outputDir = outputDir,
    maxCores = maxCores
  )
  class(externalEstimatorSettings) <- 'externalEstimatorSettings'
  return(externalEstimatorSettings)


}


#' Estimate external model performance from statistics and internal
#' test data
#'
#' @description Estimate external performance using external statistics and an internal dataset.
#' This function reweights the internal \code{z} matrix with the objective to make the weighted means as close as
#' possible to the external means represented by \code{externalStats}. Performance measures are estimated using the
#' resulting weights, the vector of actual outcomes \code{y}, and the predicted outcome probabilities \code{p}.
#'
#'
#' @param internalData a list that includes internal data and predictions with the following fields:
#'   z: a data frame of transformed feature-outcome pairs
#'   y: a vector of outcomes
#'   p: a vector of predicted outcome probabilities
#' @param externalStats a vector of means of transformed feature-outcome pairs
#'
#' Reweighing algorithm parameters:
#' @param externalEstimatorSettings an object of class \code{externalEstimatorSettings}
#' @param createEstimationLogger create a logger in outputDirectory
#'
#' @return an object of \code{estimatedExternalPerformanceFromStatistics} with the following fields:
#'
#' \itemize{
#'   \item{\code{status}:} {execution status, either \code{'Success'} or \code{'Failure'}}
#'   \item{\code{preDiagnosis}:} {Pre-diagnosis results object}
#'   \item{\code{estimationTime}:} {execution time}
#'   \item{\code{estimation}:} {a data-frame with a single column \code{'value'} and row-names that determine various
#'   estimators.}
#'   }
#'
#' @export
estimateExternalPerformanceFromStatistics <- function(
    internalData,
    externalStats,
    externalEstimatorSettings,
    createEstimationLogger = T
    )
{

  result <- list(
    status = NULL,
    preDiagnosis = NULL,
    estimationTime = NULL,
    estimation = NULL
  )
  class(result) <- 'estimatedExternalPerformanceFromStatistics'
  # Create logger
  if (createEstimationLogger) {
    analysisId <- gsub(':', '-', Sys.time())
    logFileName = file.path(externalEstimatorSettings$outputDir, glue('evaluationWithStats{analysisId}.txt'))
    logger <- ParallelLogger::createLogger(
      name = "PARALLEL",
      threshold = "TRACE",
      appenders = list(createFileAppender(
        layout = ParallelLogger::layoutParallel,
        fileName = logFileName)))
    ParallelLogger::registerLogger(logger)
    on.exit(ParallelLogger::unregisterLogger(logger))  # TODO check
  }
  # Pre diagnostics
  preD <- preDiagnostics(
    internalData$z,
    externalStats,
    maxProp = externalEstimatorSettings$maxProp,
    maxDiff = externalEstimatorSettings$maxDiff
    )
  result$preDiagnosis <- preD
  if (preD$status != 'Success') {
    ParallelLogger::logError(glue('Pre-balancing diagnosis status = {preD$status}'))
    result$status <- 'Failure'
    return(result)
  }
  # Maintain features and samples according to pre-diagnostic evaluations
  externalStats <- externalStats[preD$representedFeatures]  # Features with non-Na entries
  internalData$z <- internalData$z[preD$zidx, preD$representedFeatures]
  internalData$y <- internalData$y[preD$zidx]
  internalData$p <- internalData$p[preD$zidx]
  # Generate indices of sub-samples
  n <- sum(preD$zidx)
  nSubsets <- externalEstimatorSettings$nRepetitions
  if (externalEstimatorSettings$nMaxReweight < n*0.75)  { # TODO figure this out
    if (externalEstimatorSettings$stratified)
      idxs <- stratifiedSplitToRandomSubsets(n, externalEstimatorSettings$nMaxReweight, internalData$y, nSubsets)
    else
      idxs <- splitToRandomSubsets(n, externalEstimatorSettings$nMaxReweight, nSubsets)
  } else
    idxs <- getBootstrapSubsets(n, externalEstimatorSettings$nMaxReweight, nSubsets)

  # Estimation
  estimationStartTime <- Sys.time()
  nCores <- parallel::detectCores()
  uCores <- min(nCores-1, length(idxs), externalEstimatorSettings$maxCores)
  ParallelLogger::logInfo(glue('Detected {nCores} cores and using {uCores}\n'))
  cl <- makeCluster(uCores)
  resultsList <- clusterApply(
    cl, idxs, estimateSubsetPerformence,
    internalData=internalData, externalStats = externalStats, estimationParams=externalEstimatorSettings)
  stopCluster(cl)

  # Aggregate results into a matrix
  r1 <- resultsList[[1]]
  resultsMatrix <- matrix(nrow = nSubsets, ncol = length(r1), dimnames = list(NULL,names(r1)))
  for (k in 1:nSubsets) {
    if (!is.null(resultsList[[k]]))
      resultsMatrix[k ,] <- resultsList[[k]]
  }
  result$estimationTime <- Sys.time() - estimationStartTime
  if (externalEstimatorSettings$nMaxReweight >= n)
    meanResults <- estimateFullSetPerformance(internalData, externalStats, externalEstimatorSettings)
  else
    meanResults <- colMeans(resultsMatrix, na.rm = T)
  if (is.null(meanResults)) {
    result$status = 'Failure'
  }
  s <- summarizeBootstrap(resultsMatrix)

  # Process post diagnostic information
  # TODO add a post diagnostic that measures 'effective sample size'
  maxWSMD <- meanResults['Max Weighted SMD']
  if (is.null(maxWSMD) || is.na(maxWSMD) || maxWSMD > externalEstimatorSettings$maxWSMD) {
    ParallelLogger::logError(glue('Max weighted SMD = {maxWSMD} (Th={externalEstimatorSettings$maxWSMD})'))  # TODO which feature?
    result$status = 'Failure-large-WSMD'
  }
  else {
    wsmdWarnTh <- 0.05
    if (maxWSMD > wsmdWarnTh)
      ParallelLogger::logWarn(glue('Max weighted SMD > {wsmdWarnTh}'))  # TODO which feature?
    result$status = 'Success'
  }
  allResults <- unlist(c(list(as.list(meanResults), as.list(s))))
  result$estimation = data.frame(value=allResults)
  return(result)
}



splitToRandomSubsets <- function(n, nmax, nSubsets=1) {
  # boundaries <- getBoundaries(n, nmax)
  # nSubsets <- length(boundaries$starts)
  idxs <- vector(mode = 'list', length = nSubsets)
  nmax <- min(n, nmax)
  for (k in 1:nSubsets) {
    idxs[[k]] <- sample(n, nmax)
  }
  return(idxs)
}


getBootstrapSubsets <- function(n, nmax, nSubsets=1) {
  idxs <- vector(mode = 'list', length = nSubsets)
  nmax <- min(n, nmax)
  for (k in 1:nSubsets) {
    idxs[[k]] <- sample(n, nmax, replace = T)
  }
  return(idxs)
}



stratifiedSizes <- function(n, nmax, Y) {
  n1 = sum(Y==1)
  n0 = sum(Y==0)

  if (n1<nmax/2) {
    nmax1 <- n1
    nmax0 <- min(nmax-n1, n0)
  }
  else {
    if (n0<nmax/2) {
      nmax0 <- n0
      nmax1 <- min(nmax-n0, n1)
    }
    else {
      nmax0 <- round(nmax/2)
      nmax1 <- nmax - nmax0
    }
  }
  return(c(nmax0, nmax1))
}


stratifiedSplitToRandomSubsets <- function(n, nmax, Y, nSubsets=1) {
  n0 = sum(Y==0)
  n1 = sum(Y==1)
  nmaxs <- stratifiedSizes(n, nmax, Y)
  idxs0 <- which(Y==0)
  idxs1 <- which(Y==1)
  cat('nmax 0:', nmaxs[1], ', namx 1:', nmaxs[2], '\n')
  idxs <- vector(mode = 'list', length = nSubsets)
  for (k in 1:nSubsets) {
    idxs[[k]] <- c(idxs0[sample(n0, nmaxs[1])], idxs1[sample(n1, nmaxs[2])])
  }
  return(idxs)
}


estimateSubsetPerformence <- function(subsetIdxs, internalData, externalStats, estimationParams) {

  internalData$z <- internalData$z[subsetIdxs, ]
  internalData$y <- internalData$y[subsetIdxs]
  internalData$p <- internalData$p[subsetIdxs]

  return(estimateFullSetPerformance(internalData, externalStats, estimationParams))
}


estimateFullSetPerformance <- function(internalData, externalStats, estimationParams) {
  nZ <- nrow(internalData$z)
  pZ <- ncol(internalData$z)
  pa <- estimationParams
  ParallelLogger::logInfo(glue('Reweighting data dimensions: n={nZ}, p={pZ}'))
  ParallelLogger::logInfo(
    glue('Parameters: {pa$divergence}, lambda={pa$lambda}, distance={pa$distance}, {pa$optimizationMethod}'))

  dbRes <- list()
  y <- internalData$y
  z <- internalData$z
  p <- internalData$p

  classValues <- unique(y)
  nClasses <- length(classValues)
  if (nClasses != 2) {
    cat('Class values are', classValues, '\n')
    ParallelLogger::logError(glue('Bad namber of classes {nClasses}'))
    return (NULL)
  }
  # Re-weighting
  w <- reweightByMeans(
    z, externalStats,
    divergence = pa$divergence, lambda = pa$lambda, minSd = pa$minSd, minW = pa$minW, distance = pa$distance,
    optimizationMethod = pa$optimizationMethod, solver = pa$solver,
    verbose = T)
  if (sum(is.na(w))==0) {
    widx <- w>0
    dbRes[['n']] <- sum(widx)
    if (is.factor(y)) # TODO is this the right place
      y <- as.numeric(y)-1
    dbRes[['n outcome']] <- as.numeric(t(widx) %*% y)
    # Post diagnostics
    postD <- postDiagnostics(w, z, externalStats)
    dbRes <- c(dbRes, postD)
    # Performance measures
    m <- getPerformanceMeasures(y[widx], p[widx], w[widx])
    dbRes <- c(dbRes, m)
    return (unlist(dbRes))
  } else
    return (NULL)
}




#' Pre re-weighting diagnostics
#'
#' @description compare external expectations and internal means before reweighting
#'
#' @param z a data frame of transformed feature-outcome pairs
#' @param mu a vector of means of transformed feature-outcome pairs
#' @param maxProp maximum proportion between internal and external means to determine imbalance
#' @param maxDiff maximum difference from which to check max prop
#' @param npMinRatio minimum ratio between number of features and number of examples
#'
#' @return a named list with the following fields:
#' ...
#'
#' @export
preDiagnostics <- function(z, mu, maxProp, maxDiff, npMinRatio = 4) {
  nInput <- nrow(z)
  pInput <- ncol(z)
  ParallelLogger::logInfo(glue('Input data size n={nInput}, p={pInput}'))

  naInZ <- any(is.na(z))
  if (naInZ) {
    ParallelLogger::logError("Data set has na entries, cannot reweight")
    return(list(status='Failure'))  # Do we need this
  }
  # remove features with Na entries in table1
  naInMu <- is.na(mu)
  if (any(naInMu)) {
    ParallelLogger::logError("Expectation vector has na entries, cannot reweight")
    return(list(status='Failure'))  # Do we need this
  }
  includeFeatures <- !naInMu
  # Regardless of Na status in mu, perform other tests
  binaryResults <- preCheckBinaryFeatures(z, mu, includeFeatures, maxProp, maxDiff)
  unaryResults <- preCheckUnaryFeatures(z, mu, binaryResults, maxUnaryDiff = maxDiff)
  includeFeatures <- unaryResults$includeFeatures

  representedFeatures <- names(mu[includeFeatures])
  # Check range of variables with >1 values
  numericFeatureIndicators <- apply(z, 2, function(c) {length(unique(c))>1})
  numericFeatureIndicators <- numericFeatureIndicators & includeFeatures
  muR <- mu[numericFeatureIndicators]
  zR <- z[binaryResults$zidx, numericFeatureIndicators]  # TODO - should we limit to zidx?
  inRange <- (muR >= apply(zR, 2, min)) & (muR <= apply(zR, 2, max))
  outOfRange <- names(muR[!inRange])
  if (length(outOfRange) > 0) {
    ParallelLogger::logError('Out of range variables:')
    for (f in outOfRange)
      ParallelLogger::logError(glue('{f}, mu={muR[f]}, min={min(zR[,f])}, max={max(zR[ ,f])}'))
  }
  # few samples
  fewSamples = sum(binaryResults$zidx)/length(representedFeatures) < npMinRatio
  if (fewSamples) {
    cat('Few samples\n')
    ParallelLogger::logError(glue("Few samples n={sum(binaryResults$zidx)}, p={length(representedFeatures)}"))
  }

  status = 'Success'
  if ( length(outOfRange) > 0 || unaryResults$incompatableUnaryVariable || fewSamples
       # || length(binaryResults$highlySkewedBinary)>0
       # || length(unaryResults$unaryFeatures) > 0
  )
    status = 'Failure'
  ParallelLogger::logInfo(glue('Pre-evaluation diagnosis status = {status}'))
  return (list(
    outOfRange = outOfRange,
    representedFeatures=representedFeatures,
    zidx = binaryResults$zidx,
    unaryFeatures = unaryResults$unaryFeatures,
    incompatableUnaryVariable = unaryResults$incompatableUnaryVariable,
    highlySkewedBinary=binaryResults$highlySkewedBinary,
    status = status
    ))
}


preCheckBinaryFeatures <- function(z, mu, includeFeatures, maxProp, maxDiff) {
  n1 <- nrow(z)
  binaryFeatureIndicators <- apply(z, 2, function(c) {length(unique(c))==2})
  binaryFeatureIndicators <- binaryFeatureIndicators & includeFeatures
  # Identify feature that has a single value in the external dataset and the corresponding sub-population in the
  # internal one
  binaryFeatures <- names(mu[binaryFeatureIndicators])
  zidx <- rep(TRUE, n1)
  for (f in binaryFeatures) {
    vals <- unique(z[,f])
    for (val in vals) {
      if (mu[f]==val) {
        ParallelLogger::logInfo(glue('Removing subjects in which binary {f} != {val} and mu == {val}'))
        zidx <- zidx & (z[ , f]==val)
        includeFeatures[f] <- F
        ParallelLogger::logInfo(glue('Maintained {sum(zidx)}/{n1} subjects and {sum(includeFeatures)} features.\n'))
      }
    }
  }
  includeBinaryFeatures <- binaryFeatureIndicators & includeFeatures
  binaryFeatures <- names(mu[includeBinaryFeatures])
  ParallelLogger::logInfo(glue('z has {length(binaryFeatures)} binary features'))
  highlySkewedBinary <- getHighlySkewedBinaryFeatures(mu[binaryFeatures], z[zidx, binaryFeatures], maxProp, maxDiff)
  return(list(includeFeatures=includeFeatures, zidx = zidx, highlySkewedBinary=highlySkewedBinary))
}


# TODO compare with the version used in the Stroke use-case
#' Get highly skewed binary features
#'
#' Assuming that if the max proportion is violated than the expectations may be small and therefore checking also
#' maxDiff
#'
#' @param mu expectations
#' @param z data
#' @param maxProp max proportion
#' @param maxDiff max difference from with to check max proportion
#'
getHighlySkewedBinaryFeatures <- function(mu, z, maxProp, maxDiff) {
  ParallelLogger::logInfo(glue("Checking skewness with max proportion = {maxProp}, conditional max diff {maxDiff}"))
  imbalanced <- rep(F, length(mu))
  propM <- rep(1, length(mu))
  names(propM) <- names(mu)

  if (!is.vector(z)) {
    meanz <- colMeans(z)
    for (i in 1:(length(mu))) {
      minzi <- min(z[ , i])
      propM[i] <- abs((mu[i]-minzi)/(meanz[i]-minzi))
      imbalanced[i] <- (propM[i]>maxProp | (1/maxProp)>propM[i]) & (abs(mu[i]-meanz[i]) > maxDiff)
    }
  }
  else {
    meanz <- mean(z)
    minz <- min(z)
    propM[1] <- abs((mu[1]-minz)/(meanz-minz))
    imbalanced[1] <- (propM[1]>maxProp | (1/maxProp)>propM[1]) & (abs(mu[1]-meanz) > maxDiff)
  }
  skewedNames <- names(mu[imbalanced])
  if (sum(imbalanced) > 0) {
    reportDf <- data.frame(matrix(ncol=3, nrow=sum(imbalanced)))
    rownames(reportDf) <- names(mu[imbalanced])
    colnames(reportDf) <- c('Int', 'Ext', 'Prop')
    reportDf[[1]] <- meanz[imbalanced]
    reportDf[[2]] <- mu[imbalanced]
    reportDf[[3]] <- propM[imbalanced]
    cat('Found imbalanced features:\n')
    print(reportDf)
    for (f in skewedNames) {
      ParallelLogger::logWarn(glue('Skewed feature {f}: mean={meanz[f]} mu={mu[f]} r={propM[f]}'))
    }
  }
  return(skewedNames)
}



preCheckUnaryFeatures <- function(z, mu, results, maxUnaryDiff) {
  n1 <- nrow(z)
  includeFeatures <- results$includeFeatures
  unaryFeatures <- apply(z, 2, function(c) {length(unique(c))==1})
  unaryFeatures <- unaryFeatures & includeFeatures
  featureNames <- names(mu)

  incompatableUnaryVariable <- F
  nUnary <-sum(unaryFeatures)
  if (nUnary>0) {
    for (f in (featureNames[unaryFeatures]))
      ParallelLogger::logInfo(glue("Found unary feature {f}"))
    if (nUnary>1)
      badUnary <- abs(colMeans(z[, unaryFeatures])-mu[unaryFeatures]) > maxUnaryDiff  # TODO consider sample size
    else
      badUnary <- abs(mean(z[, unaryFeatures])-mu[unaryFeatures]) > maxUnaryDiff
    if (any(badUnary)) {
        incompatableUnaryVariable <- T
        for (f in (featureNames[unaryFeatures][badUnary]))
          ParallelLogger::logError(glue('Bad unary feature {f}, internal={mean(z[ ,f])}, external={mu[f]}'))
    }
    includeFeatures <- includeFeatures & !unaryFeatures
  }
  return(list(
    includeFeatures=includeFeatures,
    incompatableUnaryVariable = incompatableUnaryVariable,
    unaryFeatures = featureNames[unaryFeatures]))
}


#' Summarize bootstrap
#'
#' @description Summarize statistics of metrics obtained by bootstrapping
#'
#' @param b bootstrap results matrix columns correspond to different metrix, rows to repetitions. Columns
#' should be named by the metric.
#'
#' @return a named list with bootstrap statistics for every metric
#'
summarizeBootstrap <- function(b) {
  probs <- c(0.025, 0.5, 0.975)
  nboot <- nrow(b)
  s <- list()
  for (measure in colnames(b)) {
    r <- b[,measure]  # TODO learn how to extract vectors from matrices
    resultsQuantiles <- quantile(r, probs = probs, na.rm = TRUE)
    # for (i in 1:length(probs))
    s[[paste('Rough 95% lower', measure)]] = resultsQuantiles[[1]]
    s[[paste('Median', measure)]] = resultsQuantiles[[2]]
    s[[paste('Rough 95% upper', measure)]] = resultsQuantiles[[3]]

    s[[paste(measure, 'mean')]] = mean(r, na.rm = TRUE)
    s[[paste(measure, 'sd')]] = sd(r, na.rm = TRUE)
    s[['n repetitions']] = nboot - sum(is.na(r))
  }
  return(s)
}


#' Post reweighing diagnostics
#'
#' @description compute diagnostics of weighted samples
#'
#'
#' @param w a vector of weights
#' @param z a data frame of transformed feature-outcome pairs
#' @param mu a vector of means of transformed feature-outcome pairs
#'
#' @return a named list with the following:
#'   maxw
#' TODO   W2
#' kl
#'
postDiagnostics <- function(w, z, mu) {
  n <- length(w)
  p <- w/sum(w)
  klIdx <- p>0
  kl <- log(n) + as.numeric(t(p[klIdx]) %*% log(p[klIdx]))
  chi2ToUnif <- n*sum((p-1/n)**2)  #  = \sum(p-1/n)^2/1/n
  maxWeightedSMD <- computeMaxSMD(mu, z, p)

  diagnostics <- list(
    'Max weight'=max(w),
    'chi2 to uniform' = chi2ToUnif,
    kl = kl,
    'Max Weighted SMD' = maxWeightedSMD)  # TODO add diagnostics
  return(diagnostics)
}


#' Get performance measures
#'
#' @description get performance measure from a pair of binary outcome vector and model predictions. The observations
#' may be weighted by a weight vector
#'
#' @param y a binary outcome vector
#' @param p probabilities vector
#' @param w (optional) weight vector
#'
getPerformanceMeasures <- function(y, p, w=NULL) {
  nClasses <- length(unique(y))
  if (nClasses==2) {
    if (is.null(w))
      pAuc <- as.numeric(auc(roc(y, p, direction='<', quiet=T)))
    else
      pAuc <- WeightedAUC(WeightedROC(p, y, w))
    pLogLike <- WeightedLogLike(y, p, w)
    pBrier <- WeightedBrier(y, p, w)
    pMeanObservedRisk <- meanObservedRisk(y, w)
    pMeanPredictionRisk <- meanPredictionRisk(p, w)
  } else {
    warning(glue('Non binary outcome vector, number of classess = {nClasses}'))
    return(NULL)
  }
  return(list('AUROC' = pAuc,
              'Log likelyhood' = pLogLike,
              'Brier score' = pBrier,
              'Global calibration mean prediction' = pMeanPredictionRisk,
              'Global calibration observed risk' = pMeanObservedRisk
              )
         )
}

#' Estimate internal performance
#'
#' @description estimate internal performance of a model with bootstrap based confidence interval
#'
#' @param y a binary outcome vector
#' @param p probabilities vector
#' @param nboot number of bootstrap repetitions
#'
#' @export
estimateInternalPerformance <- function(y, p, nboot) {
  dbRes <- list()
  dbRes[['n']] <- length(y)
  dbRes[['n outcome']] <- sum(y)
  m <- getPerformanceMeasures(y, p)
  dbRes <- c(dbRes, m)
  br <- matrix(nrow = nboot, ncol = length(m), dimnames = list(NULL, names(m)))
  for (ib in 1:nboot) {
    ii <- sample(length(y), length(y), replace = TRUE)
    m <- getPerformanceMeasures(y[ii], p[ii])
    br[ib, ] <- unlist(m)
  }
  s <- summarizeBootstrap(br)
  dbRes <- c(dbRes, s)
  return(dbRes)
}
