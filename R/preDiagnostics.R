
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
  mu <- mu[muIntersection]
  z <- z[, muIntersection]
  if (length(muIntersection)<2)  # Y and an additional feature
    return(list(status='Failure'))

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
  binaryResults <- preCheckBinaryFeatures(z, mu, includeFeatures)
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
      ParallelLogger::logWarn(glue('{f}, mu={muR[f]}, min={min(zR[,f])}, max={max(zR[ ,f])}'))
  }
  # few samples
  fewSamples = sum(binaryResults$zidx)/length(representedFeatures) < npMinRatio
  if (fewSamples) {
    ParallelLogger::logWarn(glue("Few samples n={sum(binaryResults$zidx)}, p={length(representedFeatures)}"))
  }

  status = 'Success'
  if ( length(outOfRange) > 0 || unaryResults$incompatableUnaryVariable || fewSamples
       || length(binaryResults$highlySkewedBinary)>0
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


preCheckBinaryFeatures <- function(z, mu, includeFeatures) {
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
  highlySkewedBinary <- getHighlySkewedBinaryFeatures(mu[binaryFeatures], z[zidx, binaryFeatures])
  includeFeatures[binaryFeatures] <- highlySkewedBinary$includeFeature
  return(list(
    includeFeatures=includeFeatures,
    zidx = zidx,
    highlySkewedBinary=highlySkewedBinary$skewedNames))
}


varProxy <- function(n, n1, muExt) {
  n0 <- n-n1
  return( ((1-muExt)**2)/n0 + (muExt**2)/n1 )
}


#' Get highly skewed binary features
#'
#' @param mu expectations
#' @param z data
#' @param minNumReport minimum number of samples to report in the logger
#'
#' TODO better treatment cases in which binary values may not be 0 or 1. Transform to this form if needed
#'
getHighlySkewedBinaryFeatures <- function(mu, z, minNumReport=20, maxDiff=0.01) {
  imbalanced <- rep(F, length(mu))
  includeFeatures <- rep(T, length(mu))
  n1s <- rep(NA, length(mu))
  varProxies <- rep(NA, length(mu))

  if (is.vector(z))
    z <- as.matrix(z)
  n <- nrow(z)
  meanz <- colMeans(z)
  for (i in 1:(length(mu))) {
    minzi <- min(z[ , i])
    maxzi <- max(z[ , i])
    if (maxzi > minzi) {  # feature i is still binary after removal of samples
      n1s[i] <- sum(z[ , i] == maxzi)
      varProxies[i] <- varProxy(n, n1s[i], (mu[i]-minzi)/(maxzi-minzi))
      imbalanced[i] <- varProxies[i] > max(5/n, 0.0001)
    }
    else {
      if (minzi==0)
        n1s[i] <- 0
      else
        n1s[i] <- n
      imbalanced[i] <- abs(mu[i]-minzi) > maxDiff
      varProxies[i] <- NA
    }
  }
  if (sum(imbalanced) > 0) {
    if (F) {  # TODO Aggregated skewed features, consider this
      reportDf <- data.frame(matrix(ncol=2, nrow=sum(imbalanced)))
      rownames(reportDf) <- names(mu[imbalanced])
      colnames(reportDf) <- c('Int', 'Ext')
      reportDf[[1]] <- meanz[imbalanced]
      reportDf[[2]] <- mu[imbalanced]
    }
    for (i in which(imbalanced)) {
      f <- names(mu)[i]
      smu <- format(mu[i], digits=3)
      if (n1s[i] <= minNumReport) {
        n1sStr <- glue('<={minNumReport}')
        meanStr <- glue('<={minNumReport/n}')
      } else {
        n1sStr <- glue('={n1s[i]}')
        meanStr <- glue('={format(meanz[i], digits=3)}')
      }
      if (mu[i] < maxDiff)  # in this case both mu and n1 are small
        includeFeatures[i] <- F
      else
        ParallelLogger::logWarn(glue('Skewed feature {f}: n={n} n1{n1sStr} mean{meanStr} mu={smu}'))
    }
  }
  skewedNames <- names(mu[imbalanced & includeFeatures])
  return(list(skewedNames=skewedNames, includeFeatures=includeFeatures))
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
    if (nUnary>1)
      maxUnarySkew <- max(abs(colMeans(z[, unaryFeatures])-mu[unaryFeatures]))
    else
      maxUnarySkew <- mean(z[, unaryFeatures])-mu[unaryFeatures]
    if (maxUnarySkew > maxUnaryDiff) {
      ParallelLogger::logWarn('Incopatible unary variables:')
      incompatableUnaryVariable <- T
      for (f in (featureNames[unaryFeatures]))
        if (mean(z[ ,f]) != mu[f])
          ParallelLogger::logError(glue('unary {f}, internal={mean(z[ ,f])}, external={mu[f]}'))
    }
    includeFeatures <- includeFeatures & !unaryFeatures
  }
  return(list(
    includeFeatures=includeFeatures,
    incompatableUnaryVariable = incompatableUnaryVariable,
    unaryFeatures = featureNames[unaryFeatures]))
}


checkVariableNamesOverlap <- function(z, mu) {
  zMinusMu <- !colnames(z) %in% names(mu)
  n <- sum(zMinusMu)
  if (n>0) {
    ParallelLogger::logWarn(glue('{n} variables are in z but not in mu'))
    missingVars <- colnames(z)[zMinusMu]
    for (i in 1:n)
      ParallelLogger::logWarn(glue('{missingVars[i]} is in z but not in mu'))
  }
  muMinusZ <- !names(mu) %in% colnames(z)
  n <- sum(muMinusZ)
  if (n>0) {
    ParallelLogger::logInfo(glue('{n} variables are in mu but not in z'))  # info because we assume this case may be ok
    missingVars <- names(mu)[muMinusZ]
    for (i in 1:n)
      ParallelLogger::logInfo(glue('{missingVars[i]} is in mu but not in z'))
  }
  muIntersection <- names(mu)[!muMinusZ]
  if (length(muIntersection)<2)
    ParallelLogger::logWarn(glue('z and mu have only {length(muIntersection)} overlapping variable names'))
  return(muIntersection)
}
