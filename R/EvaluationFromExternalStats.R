#' @importFrom stats sd quantile
#' @importFrom pROC auc roc
#' @importFrom WeightedROC WeightedAUC WeightedROC
#' @import ParallelLogger
NULL


#' evaluation from external stats

#' @title Create settings for \code{estimateExternalPerformanceFromStatistics}
#'
#' @param reweightAlgorithm algorithm that attempts to reweight an internal sample to get similar expectations from an
#' external one
#' @param nMaxReweight maximum number of samples for re-weighting. If this number is smaller than the number of
#' observations, then the estimator acts on sub-samples of the data-set.
#' @param nRepetitions number of repetitions with sub-samples.
#' @param stratifiedSampling attempt to sample a balanced number of samples for both outcome groups.
#' @param maxDiff maximum difference between  external and internal fequencies in case of unary variables or
#' binary variables with high proportion.
#' @param maxWSMD maximum allowed weighted standardized mean difference
#' @param outputDir output directory for logging
#' @param maxCores maximum number of cores for parallel processing of multiple sub-samples.
#' @param shortName a short name for display puposes
#' @param warmStartAlgorithm algorithm for intialization of weights before bootstrapping
#'
#' @return
#' An object of class \code{externalEstimatorSettings}
#'
#' TODO update this
#'
#' @export
createExternalEstimatorSettings <- function(
    reweightAlgorithm,
    nMaxReweight = 20000,
    nRepetitions = 1,
    stratifiedSampling = T,
    maxDiff = 0.01,
    maxWSMD = 0.2,
    outputDir = getwd(),
    maxCores = 1,
    shortName = NULL,
    warmStartAlgorithm = NULL
)
{
  if (is.null(shortName))
    shortName <- glue('{reweightAlgorithm$shortName} n max {nMaxReweight}')
  externalEstimatorSettings <- list(
    reweightAlgorithm = reweightAlgorithm,
    nMaxReweight = nMaxReweight,
    nRepetitions = nRepetitions,
    stratifiedSampling = stratifiedSampling,
    maxDiff = maxDiff,
    maxWSMD = maxWSMD,
    outputDir = outputDir,
    maxCores = maxCores,
    shortName = shortName,
    warmStartAlgorithm = warmStartAlgorithm
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
#' @param createEstimationLogger boolearn, create a logger in outputDirectory
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
  # Initialize sub-samples
  n <- sum(preD$zidx)
  nSubsets <- externalEstimatorSettings$nRepetitions
  if (externalEstimatorSettings$nMaxReweight < n*0.75)  { # TODO figure this out
    if (externalEstimatorSettings$stratified)
      subsetSampler <- stratifiedSampler(n, externalEstimatorSettings$nMaxReweight, internalData$y)
    else
      subsetSampler <- randomSubsetSampler(n, externalEstimatorSettings$nMaxReweight, replace = F)
  } else
    subsetSampler <- randomSubsetSampler(n, externalEstimatorSettings$nMaxReweight, replace = T)  # Bootstrap

  # Estimation
  estimationStartTime <- Sys.time()
  if (!is.null(externalEstimatorSettings$warmStartAlgorithm)) {    # Warm start if defined
    externalEstimatorSettings <- warmStart(externalEstimatorSettings, internalData, externalStats)
  }
  nCores <- parallel::detectCores()
  uCores <- min(nCores-1, nSubsets, externalEstimatorSettings$maxCores)
  ParallelLogger::logInfo(glue('Detected {nCores} cores and using {uCores}\n'))
  cl <- makeCluster(uCores)
  resultsList <- clusterApply(
    cl, 1:nSubsets, estimateSubsetPerformence,
    internalData=internalData, externalStats = externalStats, estimationParams=externalEstimatorSettings, subsetSampler)
  stopCluster(cl)
  result$estimationTime <- difftime(Sys.time(), estimationStartTime, units='mins')

  resultsMatrix <- transformResultListToMatrix(resultsList)
  if (is.null(resultsMatrix)) {
    ParallelLogger::logError('All results are NULL')
    result$status <- 'Failure'
    return(result)
  }
  meanResults <- colMeans(resultsMatrix, na.rm = T)
  if (is.null(meanResults)) {
    result$status = 'Failure'
    return(result)
  }
  s <- summarizeBootstrap(resultsMatrix)

  # Process post diagnostic information
  # TODO add a post diagnostic that measures 'effective sample size'
  result$status <- checkWsmdStatus(resultsMatrix, meanResults, externalEstimatorSettings)
  result$estimation = data.frame(value=unlist(c(list(as.list(meanResults), as.list(s)))))
  return(result)
}


#' Initialize parameters using a short run on the entire set
#'
#' @param externalEstimatorSettings an object of class \code{externalEstimatorSettings} with the parameters
#' of the reweighing algorithm
#' @param internalData a list that includes internal data and predictions with the following fields:
#'   z: a data frame of transformed feature-outcome pairs
#'   y: a vector of outcomes
#'   p: a vector of predicted outcome probabilities
#' @param externalStats a vector of means of transformed feature-outcome pairs
#'
#' @return an object of class \code{externalEstimatorSettings}
#'
warmStart <- function(externalEstimatorSettings, internalData, externalStats) {
  warmSettings <- externalEstimatorSettings
  warmSettings$reweightAlgorithm <- externalEstimatorSettings$warmStartAlgorithm
  mainResult <- estimateFullSetPerformance(internalData, externalStats, warmSettings)
  w <- mainResult$reweightResults$w
  if (!is.null(mainResult)) {
    ra <- externalEstimatorSettings$reweightAlgorithm
    externalEstimatorSettings$reweightAlgorithm <- ra$setInitialValue(ra, mainResult$reweightResults)
    wRA <- externalEstimatorSettings$reweightAlgorithm$w0  # TODO just for debugging
    ParallelLogger::logInfo(glue('Reweight algorithm w0: sum={sum(wRA)}, sd={sd(wRA)} n0-={sum(wRA<=0)}'))
  }
  else
    ParallelLogger::logWarn('Failed to estimate initial value')
  return(externalEstimatorSettings)
}


transformResultListToMatrix <- function(resultsList) {
  resultsNames <- NULL
  nResults <- 0
  nSubsets <- length(resultsList)
  for (k in 1:nSubsets) {
    if (!is.null(resultsList[[k]])) {
      resultsNames <- names(resultsList[[k]])
      nResults <- nResults + 1
    }
  }
  if (nResults == 0) {
    return(NULL)
  }
  resultsMatrix <- matrix(nrow = nResults, ncol = length(resultsNames), dimnames = list(NULL, resultsNames))
  iResult <- 0
  for (k in 1:nSubsets) {
    if (!is.null(resultsList[[k]])) {
      iResult = iResult + 1
      resultsMatrix[iResult ,] <- resultsList[[k]]
    }
  }
  return(resultsMatrix)
}


checkWsmdStatus <- function(resultsMatrix, meanResults, externalEstimatorSettings) {
  maxWSMD <- meanResults['Max Weighted SMD']
  if (is.null(maxWSMD) || is.na(maxWSMD)) {
    ParallelLogger::logError('Max weighted SMD is NULL or Na')
    status = 'NULL-WSMD'
  } else {
    for (i in 1:nrow(resultsMatrix))
      ParallelLogger::logInfo(glue("Max weighted SMD {i} {resultsMatrix[i, 'Max Weighted SMD']}"))
    if (maxWSMD > externalEstimatorSettings$maxWSMD) {
      ParallelLogger::logError(glue('Max weighted SMD = {maxWSMD} (Th={externalEstimatorSettings$maxWSMD})'))
      status = 'Large-WSMD'
    }
    else {
      wsmdWarnTh <- 0.05
      if (maxWSMD > wsmdWarnTh)
        ParallelLogger::logWarn(glue('Max weighted SMD = {maxWSMD} > {wsmdWarnTh}'))  # TODO which feature?
      status = 'Success'
    }
  }
  return(status)
}


######################################################################################################################
#
# Subsets samplers

#' Random subset sampler
#'
#' @param n number of samples
#' @param nmax maximum number of subsamples
#' @param replace boolean for sampling with replacement
#'
#' @return a \code{randomSubsetSampler} objects
#'
randomSubsetSampler <- function(n, nmax, replace=F) {
  nmax <- min(n, nmax)

  l <- list(
    n = n,
    nmax = nmax,
    replace = replace,
    sampleNext = randomSubset
  )
  class(l) <- 'randomSubsetSampler'
  return(l)

}

randomSubset <- function(s) {
  return(sample(s$n, s$nmax, replace=s$replace))
}


#' Stratified subset sampler
#'
#' @param n number of samples
#' @param nmax maximum number of samples
#' @param Y outcome vector
#'
#' @return a \code{stratifiedSampler} objects
#'
stratifiedSampler <- function(n, nmax, Y) {
  l <- list(
    n0 = sum(Y==0),
    n1 = sum(Y==1),
    nmaxs = stratifiedSizes(n, nmax, Y),
    idxs0 = which(Y==0),
    idxs1 = which(Y==1),
    sampleNext = stratifiedSplitToRandomSubsets
  )
  class(l) <- 'stratifiedSampler'
  ParallelLogger::logInfo(glue('n0={l$n0}, n1={l$n1}, nmaxs1={l$nmaxs[1]}, nmaxs2={l$nmaxs[2]}'))
  return(l)

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


stratifiedSplitToRandomSubsets <- function(s) {
  # Replace is TRUE by default because the smaller set may not suffice to give different subset
  ParallelLogger::logInfo(glue('n0={s$n0}, n1={s$n1}, nmaxs1={s$nmaxs[1]}, nmaxs2={s$nmaxs[2]}'))
  idxs <- c(s$idxs0[sample(s$n0, s$nmaxs[1], replace = T)], s$idxs1[sample(s$n1, s$nmaxs[2], replace = T)]) # TODO examine
  return(idxs)
}


estimateSubsetPerformence <- function(i, internalData, externalStats, estimationParams, subsetSampler) {

  subsetIdxs <- subsetSampler$sampleNext(subsetSampler)
  internalData$z <- internalData$z[subsetIdxs, ]
  internalData$y <- internalData$y[subsetIdxs]
  internalData$p <- internalData$p[subsetIdxs]

  if (!is.null(estimationParams$reweightAlgorithm$w0))
    estimationParams$reweightAlgorithm$w0 <- estimationParams$reweightAlgorithm$w0[subsetIdxs]

  r <- estimateFullSetPerformance(internalData, externalStats, estimationParams)

  if (!is.null(r))
    return(r$dbRes)
  else
    return(NULL)  # TODO check what happens in this case
}


estimateFullSetPerformance <- function(internalData, externalStats, estimationParams) {
  nZ <- nrow(internalData$z)
  pZ <- ncol(internalData$z)
  pa <- estimationParams
  ParallelLogger::logInfo(glue('Reweighting data dimensions: n={nZ}, p={pZ}'))
  # ParallelLogger::logInfo(
  #   glue('Parameters: {pa$divergence}, lambda={pa$lambda}, distance={pa$distance}, {pa$optimizationMethod}'))

  dbRes <- list()

  classValues <- unique(internalData$y)
  nClasses <- length(classValues)
  if (nClasses != 2) {
    ParallelLogger::logError(glue('Bad namber of classes {nClasses}, values {classValues}'))
    return (NULL)
  }
  # Re-weighting
  wOptimizer <- pa$reweightAlgorithm
  reweightResults <- wOptimizer$optimize(wOptimizer, internalData$z, externalStats)
  dbRes[['Opt err']] <- reweightResults$err
  dbRes[['n iter']] <- nrow(reweightResults$log)

  # TODO make this optional
  if (F) {
    runId <- gsub(':', '-', Sys.time())
    logName = glue('{class(wOptimizer)} {runId} log.png')
    if (class(wOptimizer) %in% c('sgdWeightOptimizer', 'sgdTunedWeightOptimizer', 'seTunedWeightOptimizer')) {
      png(filename = file.path(wOptimizer$outputDir, logName), width = 4*480, height = 4*480)
      plotOptimizationLog(reweightResults$log, wOptimizer, nrow(internalData$z))
      dev.off()
    }
  }

  w <- reweightResults$w
  if (sum(is.na(w))==0) {
    widx <- w>0
    dbRes[['n']] <- sum(widx)
    if (is.factor(internalData$y)) # TODO is this the right place
      internalData$y <- as.numeric(internalData$y)-1
    dbRes[['n outcome']] <- as.numeric(t(widx) %*% internalData$y)
    # Post diagnostics
    postD <- postDiagnostics(w, internalData$z, externalStats)
    dbRes <- c(dbRes, postD)
    # Performance measures
    m <- getPerformanceMeasures(internalData$y[widx], internalData$p[widx], w[widx])
    dbRes <- c(dbRes, m)
    return (list(dbRes=(unlist(dbRes)), reweightResults=reweightResults))
  } else {
    ParallelLogger::logWarn('w contains NA')
    return (NULL)
  }
}


checkVariableNamesOverlap <- function(z, mu) {
  zMinusMu <- !colnames(z) %in% names(mu)
  n <- sum(zMinusMu)
  if (n>0) {
    ParallelLogger::logInfo(glue('{n} variables are in z but not in mu'))
    missingVars <- colnames(z)[zMinusMu]
    for (i in 1:n)
      ParallelLogger::logInfo(glue('{missingVars[i]} is in z but not in mu'))
  }
  muMinusZ <- !names(mu) %in% colnames(z)
  n <- sum(muMinusZ)
  if (n>0) {
    ParallelLogger::logInfo(glue('{n} variables are in mu but not in z'))
    missingVars <- names(mu)[muMinusZ]
    for (i in 1:n)
      ParallelLogger::logInfo(glue('{missingVars[i]} is in mu but not in z'))
  }
  muIntersection <- !muMinusZ
  if (sum(muIntersection)<2)
    ParallelLogger::logFatal(glue('z and mu have only {sum(muIntersection)} overlapping variable names'))
  return(muIntersection)
}



#' Pre re-weighting diagnostics
#'
#' @description compare external expectations and internal means before reweighting
#'
#' @param z a data frame of transformed feature-outcome pairs
#' @param mu a vector of means of transformed feature-outcome pairs
#' @param maxDiff maximum difference from which to check max prop
#' @param npMinRatio minimum ratio between number of features and number of examples
#'
#' @return a named list with the following fields:
#' ...
#'
#' @export
preDiagnostics <- function(z, mu, maxDiff, npMinRatio = 4) {
  nInput <- nrow(z)
  pInput <- ncol(z)

  muIntersection <- checkVariableNamesOverlap(z, mu)

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
  includeFeatures <- (!naInMu) & muIntersection
  # Regardless of Na status in mu, perform other tests
  binaryResults <- preCheckBinaryFeatures(z, mu, includeFeatures, maxDiff)
  unaryResults <- preCheckUnaryFeatures(z, mu, binaryResults, maxUnaryDiff = maxDiff)
  includeFeatures <- unaryResults$includeFeatures

  representedFeatures <- names(mu[includeFeatures])
  # Check range of variables with >1 values
  numericFeatureIndicators <- apply(z, 2, function(c) {length(unique(c))>1})
  numericFeatureIndicators <- numericFeatureIndicators & includeFeatures
  muR <- mu[numericFeatureIndicators]
  zR <- z[binaryResults$zidx, numericFeatureIndicators]  # TODO - should we limit to zidx?
  inRangeTol <- maxDiff  # TODO should we use a different param?
  inRange <- (muR >= apply(zR, 2, min)-maxDiff) & (muR <= apply(zR, 2, max)+maxDiff)
  outOfRange <- names(muR)[!inRange]
  if (length(outOfRange) > 0) {
    ParallelLogger::logError('Out of range variables:')
    for (f in outOfRange)
      ParallelLogger::logError(glue('{f}, mu={muR[f]}, min={min(zR[,f])}, max={max(zR[ ,f])}'))
  }
  # few samples
  fewSamples = sum(binaryResults$zidx)/length(representedFeatures) < npMinRatio
  if (fewSamples) {
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


preCheckBinaryFeatures <- function(z, mu, includeFeatures, maxDiff) {
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
  highlySkewedBinary <- getHighlySkewedBinaryFeatures(mu[binaryFeatures], z[zidx, binaryFeatures], maxDiff)
  return(list(includeFeatures=includeFeatures, zidx = zidx, highlySkewedBinary=highlySkewedBinary))
}


varProxy <- function(n, n1, muExt) {
  n0 <- n-n1
  return( ((1-muExt)**2)/n0 + (muExt**2)/n1 )
}


# TODO compare with the version used in the Stroke use-case
#' Get highly skewed binary features
#'
#' Assuming that if the max proportion is violated than the expectations may be small and therefore checking also
#' maxDiff
#'
#' @param mu expectations
#' @param z data
#' @param maxDiff max difference from with to check max proportion
#'
getHighlySkewedBinaryFeatures <- function(mu, z, maxDiff) {
  imbalanced <- rep(F, length(mu))
  propM <- rep(NA, length(mu))
  n1s <- rep(NA, length(mu))
  varProxies <- rep(NA, length(mu))

  if (!is.vector(z)) {
    n <- nrow(z)
    meanz <- colMeans(z)
    for (i in 1:(length(mu))) {
      minzi <- min(z[ , i])
      maxzi <- max(z[ , i])
      n1s[i] <- sum(z[ , i] == maxzi)
      propM[i] <- abs((mu[i]-minzi)/(meanz[i]-minzi))
      varProxies[i] <- varProxy(n, n1s[i], (mu[i]-minzi)/(maxzi-minzi))
      imbalanced[i] <- varProxies[i] > 0.001
    }
  }
  else {
    n <- length(z)
    meanz <- mean(z)
    minz <- min(z)
    maxz <- max(z)
    n1 <- sum(z==maxz)
    propM[1] <- abs((mu[1]-minz)/(meanz-minz))
    varProxies[1] <- varProxy(n, n1, (mu[1]-minz)/(maxz-minz))
    imbalanced[1] <- varProxies[1] > 0.001
  }
  skewedNames <- names(mu[imbalanced])
  if (sum(imbalanced) > 0) {
    reportDf <- data.frame(matrix(ncol=3, nrow=sum(imbalanced)))
    rownames(reportDf) <- names(mu[imbalanced])
    colnames(reportDf) <- c('Int', 'Ext', 'Prop')
    reportDf[[1]] <- meanz[imbalanced]
    reportDf[[2]] <- mu[imbalanced]
    reportDf[[3]] <- propM[imbalanced]
    # print(reportDf)
    for (i in which(imbalanced)) {
      f <- names(mu)[i]
      smean <- format(meanz[i], digits=3)
      smu <- format(mu[i], digits=3)
      spropM <- format(propM[i], digits=3)
      sv <-format(varProxies[i], digits=3)
      ParallelLogger::logWarn(glue('Skewed feature {f}: n={n} n1={n1s[i]} mean={smean} mu={smu} r={spropM} var proxy={sv}'))
    }
  }
  return(skewedNames)
}



preCheckUnaryFeatures <- function(z, mu, results, maxUnaryDiff=0.01) {
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
#' see https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package/
#'
#' @param b bootstrap results matrix columns correspond to different metrix, rows to repetitions. Columns
#' should be named by the metric.
#'
#' @return a named list with bootstrap statistics for every metric
#'
#' @export
summarizeBootstrap <- function(b) {
  probs <- c(0.025, 0.975)
  nboot <- nrow(b)
  s <- list()
  for (measure in colnames(b)) {
    r <- b[,measure]  # TODO learn how to extract vectors from matrices
    if (min(r)>-1) {

      resultsQuantiles <- quantile(r, probs = probs, na.rm = TRUE)
      s[[paste('prc lower', measure)]] = resultsQuantiles[[1]]
      s[[paste('prc upper', measure)]] = resultsQuantiles[[2]]

      # Standard normal interval assuming a zero bias
      # TODO correct for biases in case there is a complete sample
      # TODO tailor transformations to specific metrics
      r1 <- log(1+r)
      m <- mean(r1)
      se <- sd(r1)

      s[[paste('95% lower', measure)]] = exp(m-1.96*se)-1
      s[[paste('Median', measure)]] = median(r)
      s[[paste('95% upper', measure)]] = exp(m+1.96*se)-1

      s[[paste(measure, 'mean')]] = mean(r, na.rm = TRUE)
      s[[paste(measure, 'sd')]] = sd(r, na.rm = TRUE)
    }
  }
  s[['n repetitions']] = nboot - sum(is.na(r))
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
#' \itemize{
#' \item {\code{Max weight}:} {Maximal weight}
#' \item {\code{chi2ToUnif}:} {Chi squared divergence between weights and a uniform distribution}
#' \item {\code{kl}:} {KL distance to uniform distribution}
#' \item {\code{Max Weighted SMD}:} {Maximum over features of SMD between internal set and external statistics}
#' }
#'
#' @export
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
  ParallelLogger::logInfo(glue("Max weighted SMD = {maxWeightedSMD}"))
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
#' @return Performence measures....
#'
#' @export
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
