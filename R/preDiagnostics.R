#' Get field names of structured pre-diagnostic object
#'
#' @return a character vector of field names
#'
#' @export
getPreDiagnosticsFieldNames <- function() {
  fields <- c(
    "preDiagnosisStatus", "missingInMu", "missingInZ", "nOverlapFeatures", "naInZ", "naInMu",
    "highlySkewedBinary", "highSkewRepresentative", "incompatableUnaryVariable", "incompatibleUnaryRepresentative",
    "outOfRange", "outOfRangeRepresentative", "fewSamples", "fewSamplesDescription"
  )
  return(fields)
}


#' Pre re-weighting diagnostics
#'
#' @description compare external expectations and internal means before reweighting
#'
#' The following tests are made:
#' 1. z and mu have the same set of variables
#' 2. z and mu do not contain NA entries
#' 3. skew in binary features
#'
#' Additionally, internal samples are removed when ...
#'
#' @param z a data frame of transformed feature-outcome pairs
#' @param mu a vector of means of transformed feature-outcome pairs
#' @param maxDiff maximum difference from which to check max prop
#' @param npMinRatio minimum ratio between number of features and number of examples
#' @param maxSubset=20000
#'
#' @return a named list with the following fields:
#' \itemize{
#'   \item{outOfRange}{list of external statistics is out of range of a numeric feature}
#'   \item{representedFeatures}{names of features }
#'   \item{zidx}{indices of samples that should be used in reweighting}
#'   \item{highlySkewedBinary}{list of highly skewed variables}
#'   \item{status}{a string that determines}
#'   \item{structuredLog}{A data frame with indicators of prediagnosis results and descriptions, see below}
#' }
#'
#' Structured log is a data-frame with a single value column and the following row names:
#' \describe{
#' \item{preDiagnosisStatus}{over all status of diagnostics}
#' \item{missingInMu}{there are variables in the internal data matrix that are missing in the external statistics}
#' \item{missingInZ}{variable in the external statistics are missing in the internal data}
#' \item{nOverlapFeatures}{number of overlapping features}
#' \item{naInZ}{NA values in the internal data}
#' \item{naInMu}{NA values in the statistics vector}
#' \item{highlySkewedBinary}{highly skewed binary variables, generally indicating a few sample in the internal
#'  dataset but a high proportion of them in the external ones}
#' \item{highSkewRepresentative'}{an example for a highly skewed variable}
#' \item{incompatableUnaryVariable'}{value of an internal unary variable is different than the external one more than
#' a threshold}
#' \item{incompatibleUnaryRepresentative'}{example for incompatible unary variable}
#' \item{outOfRange}{externa statistic is out of range}
#' \item{outOfRangeRepresentative}{example}
#' \item{fewSamples}{few sample relative to number of }
#' \item{fewSamplesDescription}{elaboration on number of samples and dimension}
#' }
#'
#' @export
preDiagnostics <- function(z, mu, maxDiff, npMinRatio = 4, maxSubset=20000) {

  # Init structured log
  fields <- getPreDiagnosticsFieldNames()
  structuredLog = data.frame(value=rep(NA, length(fields)), row.names = fields)

  nInput <- nrow(z)
  pInput <- ncol(z)

  # Verify overlapping names between z columns and mu
  overlapLog <- checkVariableNamesOverlap(z, mu)
  structuredLog[names(overlapLog), 'value'] = unlist(overlapLog)
  status <- overlapLog$preDiagnosisStatus

  if (status != 'Success') {
    ParallelLogger::logError(glue('variable name overlap status = {status}'))
    return(list(structuredLog=structuredLog, status='Failure'))
  }

  # Once we verified names overlap, we align them
  mu <- mu[colnames(z)]

  # Check that there are no Na entries in z
  structuredLog['naInZ', 'value'] <- any(is.na(z))
  if (structuredLog['naInZ', 'value']) {
    ParallelLogger::logError("Data set has na entries, cannot reweight")
    status <- 'Failure'
    structuredLog$preDiagnosisStatus <- 'Failure'
    return(list(structuredLog=structuredLog, status='Failure'))
  }
  # remove features with Na entries in table1
  structuredLog['naInMu', 'value'] <- any(is.na(mu))
  if (structuredLog['naInMu', 'value']) {
    ParallelLogger::logError("Expectation vector has na entries, cannot reweight")
    status <- 'Failure'
    structuredLog$preDiagnosisStatus <- 'Failure'
    return(list(structuredLog=structuredLog, status='Failure'))
  }

  # Regardless of Na status in mu, perform other tests
  binaryResults <- preCheckBinaryFeatures(z, mu, maxSubset = maxSubset)
  structuredLog['highlySkewedBinary', 'value'] <- length(binaryResults$highlySkewedBinary)>0
  structuredLog['highSkewRepresentative', 'value'] <- binaryResults$highSkewRepresentative
  if (structuredLog['highlySkewedBinary', 'value'])
    status <- 'Failure'

  unaryResults <- preCheckUnaryFeatures(z, mu, binaryResults, maxUnaryDiff = maxDiff)
  includeFeatures <- unaryResults$includeFeatures
  structuredLog['incompatableUnaryVariable', 'value'] <- unaryResults$incompatableUnaryVariable
  structuredLog['incompatibleUnaryRepresentative', 'value'] <- unaryResults$incompatibleRepresentative
  if (structuredLog['incompatableUnaryVariable', 'value'])
    status <- 'Failure'

  representedFeatures <- names(mu[includeFeatures])
  # Check range of variables with >1 values
  numericFeatureIndicators <- apply(z, 2, function(c) {length(unique(c))>1})
  numericFeatureIndicators <- numericFeatureIndicators & includeFeatures
  muR <- mu[numericFeatureIndicators]
  zR <- z[binaryResults$zidx, numericFeatureIndicators]  # TODO - should we limit to zidx?
  inRangeTol <- maxDiff  # TODO should we use a different param?
  inRange <- (muR >= apply(zR, 2, min)-maxDiff) & (muR <= apply(zR, 2, max)+maxDiff)
  outOfRange <- names(muR)[!inRange]
  structuredLog['outOfRange', 'value'] <- length(outOfRange) > 0
  if (structuredLog['outOfRange', 'value']) {
    ParallelLogger::logError('Out of range variables:')
    f <- outOfRange[1]
    structuredLog['outOfRangeRepresentative', 'value'] <- glue('{f}, mu={muR[f]}, min={min(zR[,f])}, max={max(zR[ ,f])}')
    for (f in outOfRange)
      ParallelLogger::logWarn(glue('{f}, mu={muR[f]}, min={min(zR[,f])}, max={max(zR[ ,f])}'))
    status <- 'Failure'
  }

  # few samples
  structuredLog['fewSamples', 'value'] <- sum(binaryResults$zidx)/length(representedFeatures) < npMinRatio
  if (structuredLog['fewSamples', 'value']) {
    structuredLog['fewSamplesDescription', 'value'] <-
      glue("Few samples n={sum(binaryResults$zidx)}, p={length(representedFeatures)}")
    ParallelLogger::logWarn(structuredLog['fewSamplesDescription', 'value'])
    status <- 'Failure'
  }

  if ( status != 'Success')
    ParallelLogger::logWarn(glue('Pre-evaluation diagnosis status = {status}'))
  else
    ParallelLogger::logInfo(glue('Pre-evaluation diagnosis status = {status}'))
  structuredLog['preDiagnosisStatus', 'value'] <- status

  return (list(
    outOfRange = outOfRange,  # External statistics is out of range of a numeric feature
    representedFeatures=representedFeatures,
    zidx = binaryResults$zidx,
    highlySkewedBinary=binaryResults$highlySkewedBinary,  # list
    status = status,  # For compatibility and ease of access
    structuredLog = structuredLog
  ))
}

#' pre check binary features
#'
#' @param z data frame of features
#' @param mu named vectors of statistics
#' @param includeFeatures indicators
#' @param maxSubset max number of samples in training subsets
#'
#' @return a list
#'
preCheckBinaryFeatures <- function(z, mu, maxSubset=20000) {
  n1 <- nrow(z)
  binaryFeatureIndicators <- apply(z, 2, function(c) {length(unique(c))==2})
  includeFeatures <- rep(T, length(binaryFeatureIndicators))
  names(includeFeatures) <- names(binaryFeatureIndicators)
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
  highlySkewedBinary <- getHighlySkewedBinaryFeatures(
    mu[binaryFeatures], z[zidx, binaryFeatures], maxSubset = maxSubset)
  includeFeatures[binaryFeatures] <- highlySkewedBinary$includeFeatures
  return(list(
    includeFeatures=includeFeatures,
    zidx = zidx,
    highlySkewedBinary=highlySkewedBinary$skewedNames,
    highSkewRepresentative=highlySkewedBinary$highSkewRepresentative
  ))
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
#' @param maxDiff max difference in the unary case
#' @param maxSubset assumption about the maximum subset size
#'
#' TODO better treatment cases in which binary values may not be 0 or 1. Transform to this form if needed
#'
#' @return a list
#'
getHighlySkewedBinaryFeatures <- function(mu, z, minNumReport=20, maxDiff=0.01, maxSubset=20000) {

  imbalanced <- rep(F, length(mu))
  includeFeatures <- rep(T, length(mu))
  n1s <- rep(NA, length(mu))
  varProxies <- rep(NA, length(mu))

  if (is.vector(z))
    z <- as.matrix(z)
  n <- nrow(z)
  f <- min(maxSubset/n, 1)  # Factor for varProxy
  meanz <- colMeans(z)
  for (i in 1:(length(mu))) {
    minzi <- min(z[ , i])
    maxzi <- max(z[ , i])
    if (maxzi > minzi) {  # feature i is still binary after removal of samples
      n1s[i] <- sum(z[ , i] == maxzi)
      varProxies[i] <- varProxy(n, n1s[i]*f, (mu[i]-minzi)/(maxzi-minzi))
      if (n>2000)
        imbalanced[i] <- varProxies[i] > 0.000625 # 1/40**2
      else
        imbalanced[i] <- varProxies[i] > 2/n
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
    skewDescriptions <- rep('', length(imbalanced))
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
      if (mu[i] < maxDiff)  { # in this case both mu and n1 are small
        includeFeatures[i] <- F
        skewDescriptions[i] <- glue('Nearly-unary {f}: n={n} n1{n1sStr} mean{meanStr} mu={smu}')
      }
      else
        skewDescriptions[i] <- glue('Skewed {f}: n={n} n1{n1sStr} mean{meanStr} mu={smu}')
      ParallelLogger::logWarn(skewDescriptions[i])
    }
    varProxies[is.na(varProxies)] <- 0
    highSkewRepresentative <- skewDescriptions[which.max(varProxies)]
  }
  else
    highSkewRepresentative <- NA
  skewedNames <- names(mu[imbalanced & includeFeatures])
  # Extract high skew representative
  return(list(skewedNames=skewedNames, includeFeatures=includeFeatures, highSkewRepresentative=highSkewRepresentative))
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
    meanUnaryMu <- mu[unaryFeatures]
    if (nUnary>1) {
      meanUnaryZ <- colMeans(z[, unaryFeatures])
      unarySkew <- abs(meanUnaryZ-meanUnaryMu)
      maxUnarySkew <- max(unarySkew)
      f <- featureNames[unaryFeatures][which.max(unarySkew)]
      representative <- glue('{f}: mean z={meanUnaryZ[f]}, mu={meanUnaryMu[f]}')
    }
    else {
      meanUnaryZ <- mean(z[, unaryFeatures])
      maxUnarySkew <- abs(meanUnaryZ-meanUnaryMu)
      representative <- glue('{names(mu)[unaryFeatures]}: mean z={meanUnaryZ}, mu={meanUnaryMu}')
    }
    cat('Max unary skew', maxUnarySkew, 'max diff', maxUnaryDiff, '\n')
    if (maxUnarySkew > maxUnaryDiff) {
      ParallelLogger::logWarn('Incopatible unary variables:')
      incompatableUnaryVariable <- T
      for (f in (featureNames[unaryFeatures]))
        if (mean(z[ ,f]) != mu[f])
          ParallelLogger::logError(glue('unary {f}, internal={mean(z[ ,f])}, external={mu[f]}'))
    } else
      representative <- NA
    includeFeatures <- includeFeatures & !unaryFeatures
  } else
    representative <- NA
  return(list(
    includeFeatures=includeFeatures,
    incompatableUnaryVariable = incompatableUnaryVariable,
    unaryFeatures = featureNames[unaryFeatures],
    incompatibleRepresentative = representative))
}

#' Check overlap of variables names
#'
#' @param z a data-frame
#' @param mu named vector
#'
checkVariableNamesOverlap <- function(z, mu) {

  zMinusMu <- !colnames(z) %in% names(mu)
  n <- sum(zMinusMu)
  result <- list(preDiagnosisStatus = 'Success', missingInMu = F, missingInZ = F, nOverlapFeatures = NaN)
  if (n>0) {
    ParallelLogger::logWarn(glue('{n} variables are in z but not in mu'))
    missingVars <- colnames(z)[zMinusMu]
    for (i in 1:n)
      ParallelLogger::logWarn(glue('{missingVars[i]} is in z but not in mu'))
    result$missingInMu <- T
    result$preDiagnosisStatus <- 'Failure'
  }
  muMinusZ <- !names(mu) %in% colnames(z)
  n <- sum(muMinusZ)
  if (n>0) {
    ParallelLogger::logWarn(glue('{n} variables are in mu but not in z'))  # info because we assume this case may be ok
    missingVars <- names(mu)[muMinusZ]
    for (i in 1:n)
      ParallelLogger::logWarn(glue('{missingVars[i]} is in mu but not in z'))
    result$missingInZ <- T
    result$preDiagnosisStatus <- 'Failure'
  }
  result$nOverlapFeatures <- sum(!muMinusZ)
  if (result$nOverlapFeatures==1)
    ParallelLogger::logWarn(
      glue('z and mu have only one overlapping variable names'))
  else {
    if (result$nOverlapFeatures==0) {
      ParallelLogger::logWarn(
        glue('z and mu do not have overlapping variable names'))
      result$preDiagnosisStatus = 'Failure'
    }
  }
  return(result)
}
