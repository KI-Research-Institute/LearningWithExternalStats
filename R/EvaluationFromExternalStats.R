#' @importFrom stats sd quantile
#' @importFrom pROC auc roc
#' @importFrom WeightedROC WeightedAUC WeightedROC
#' @import ParallelLogger
NULL


#' evaluation from external stats

# TODO - replace or requires installation of ParallelLoger?

#' Estimate external performance from statistics
#'
#' @description Estimate external performance using external statistics and an internal dataset
#'
#' @param internalData a list that includes internal data and predictions with the following fields:
#'   z: a data frame of transformed feature-outcome pairs
#'   y: a vector of outcomes
#'   p: a vector of predicted outcome probabilities
#' @param externalStats a vector of means of transformed feature-outcome pairs
#'
#' Reweighing algorithm parameters:
#' @param divergence 'entropy' or 'chi2'.
#' 'entropy' directs the algorithm to minimize the negative entropy, \eqn{-\sum_i w_i \log w_i}.
#' 'chi2' is \eqn{\sum_i{w_i-\frac{1}{n}}**2}
#' @param lambda lambda - regularization parameter
#' @param minSd minimum variance for a columns to be included
#' @param optimizationMethod primal or dual. Currently dual works only with entropy divergence
#' @param maxProp maximum proportion between external and internal means
#' @param nboot number of bootstrap repetitions for confidence interval assessment
#'
#' @return a named list with the following fields:
#'   summary: a name list with pre-weighting diagnostics, post-weighting diagnostics, estimated performance measures
#'     and bootstrap based performance measures.
#'   preDiagnosis: more detailed pre-reweighing diagnosis information.
#'   bootstrap: a data frame of results of single bootstrap repetitions
#'
#' This function reweights the internal z matrix with the objective to make the weighted means as close as
#' possible to the external means represented by mu. Performance measures are estimated using the resulting
#' weights, the vector of actual outcomes y, and the predicted outcome probabilities p.
#'
#'
estimateExternalPerformanceFromStats <- function(
    internalData, externalStats,
    divergence = "entropy", lambda = 1e-06, minSd = 1e-04,
    minW = 1e-06, distance = "l2", optimizationMethod = "primal",
    maxProp = 500, nboot=10) {
  #' check the type of internal data, if in plp format, convert to
  dbRes <- list()
  # Pre diagnostics
  preD <- preDiagnostics(internalData$z, externalStats, maxProp)
  mu <- externalStats[preD$representedFeatures]  # Features with non-Na entries
  z <- internalData$z[preD$zidx, preD$representedFeatures]
  y <- internalData$y[preD$zidx]
  p <- internalData$p[preD$zidx]
  nClasses <- length(unique(y))
  if (nClasses != 2)
    stop(glue('Bad namber of classes {nClasses}'))
  if (sum(preD$numImbalanced)>0) {
    warning('imbalanced features, cannot reweight')
    return(list(summary=NULL, preDiagnosis=preD))
  }
  # Re-weighting
  w <- reweightByMeans(
    z, mu,
    divergence = divergence, lambda = lambda, minSd = minSd, minW = minW, distance = distance,
    optimizationMethod = optimizationMethod,
    verbose = T)
  if (sum(is.na(w))==0) {
    widx <- w>0
    dbRes[['n']] <- sum(widx)
    if (is.factor(y)) # TODO is this the right place
      y <- as.numeric(y)-1
    dbRes[['n y']] <- t(widx) %*% y
    # Post diagnostics
    postD <- postDiagnostics(w, z, mu)
    dbRes <- c(dbRes, postD)
    # Performance measures
    m <- getPerformanceMeasures(y[widx], p[widx], w[widx])
    dbRes <- c(dbRes, m)
    # Bootstrap
    if (nboot>0) {
      br <- matrix(nrow = nboot, ncol = length(m), dimnames = list(NULL, names(m)))
      for (ib in 1:nboot) {
        ii <- sample(nrow(z), nrow(z), replace = TRUE)
        wi <- reweightByMeans(z[ii, ], mu, lambda = lambda, verbose = F)
        if (sum(is.na(wi))==0) {
          widx <- wi>0
          m <- getPerformanceMeasures(y[ii][widx], p[ii][widx], wi[widx])
          if (!is.null(m))
            br[ib, ] <- unlist(m)
        }
      }
      s <- summarizeBootstrap(br)
      dbRes <- c(dbRes, s)
      return(list(summary=dbRes, preDiagnosis=preD, bootstrap=br))
    }
    else
      return(list(summary=dbRes, preDiagnosis=preD))
  } else
    return(list(summary=NULL, preDiagnosis=preD))
}



#' Pre re-weighting diagnostics
#'
#' @description compare external expectations and internal means before reweighting
#'
#' @param z a data frame of transformed feature-outcome pairs
#' @param mu a vector of means of transformed feature-outcome pairs
#' @param maxProp maximum proportion between internal and external means to determine imbalance
#'
#' @return a named list with the following fields:
#'   representedfeatured: a vecotr of represented feature names
#'   zidx: indices of valid subjects
#'   imbalanced: indicators of imbalance
preDiagnostics <- function(z, mu, maxProp, verbose=T) {
  # remove features with Na entries in table1
  representedFeatures <- names(mu[(!is.na(mu))] )
  ParallelLogger::logInfo(glue("Dataset has {length(representedFeatures)} non-na entries"), "\n")
  n1 <- nrow(z)
  removeUnrepresentedSubjects = T
  # TODO the following code assumes all entries represent proportions
  if (sum(mu[representedFeatures]==0)>0 & removeUnrepresentedSubjects) {
    ParallelLogger::logInfo(
      paste('Removing subjects with unrepresented propertis:', paste(mu[mu==0], collapse = ' '), sep = ' '))
    zidx <- rowSums(abs(z[, mu==0]))==0
    representedFeatures <- representedFeatures[mu[representedFeatures] != 0]
    ParallelLogger::logInfo(glue('Maintained {sum(zidx)}/{n1} subjects.\n'))
  } else {
    cat('Number of rows', n1, '\n')
    zidx = rep(TRUE, n1)
  }

  meanz <- colMeans(z[zidx, representedFeatures])
  imbalanced <- getImbalancedFeatures(mu[representedFeatures], meanz, maxProp, verbose)
  return (list(representedFeatures=representedFeatures, zidx = zidx, imbalanced=imbalanced))
}


#' Summarize bootstrap
#'
#' @description Summarize statistics of metrics obtained by bootstrapping
#'
#' @param b bootstrap results matrix columns correspond to different metrix, rows to repetitions. Columns
#' should be named by the metric.
#' @param probs quatile probabilities
#'
#' @return a named list with bootstrap statistics for every metric
#'
summarizeBootstrap <- function(b, probs=c(0.025, 0.5, 0.975)) {
  nboot <- nrow(b)
  s <- list()
  for (measure in colnames(b)) {
    r <- b[,measure]  # TODO learn how to extract vectors from matrices
    resultsQuantiles <- quantile(r, probs = probs, na.rm = TRUE)
    for (i in 1:length(probs))
      s[[paste(measure, as.character(probs[i]))]] = resultsQuantiles[[i]]

    s[[paste(measure, 'mean')]] = mean(r, na.rm = TRUE)
    s[[paste(measure, 'sd')]] = sd(r, na.rm = TRUE)
    s[[paste(measure, 'n boot')]] = nboot - sum(is.na(r))
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
  p <- w/length(w)
  klIdx <- p>0
  kl <- log(n) + t(p[klIdx]) %*% log(p[klIdx])
  chi2ToUnif <- n*sum((p-1/n)**2)  #  = \sum(p-1/n)^2/1/n
  maxWeightedSMD <- computeMaxSMD(mu, z, p)

  diagnostics <- list(maxw=max(w), 'chi2-u' = chi2ToUnif, kl = kl, maxWeightedSMD=maxWeightedSMD)  # TODO add diagnostics
  return(diagnostics)
}



# TODO specific to stroke?
getImbalancedFeatures <- function(mu, meanz, maxProp, verbose=T) {
  propM <- abs(mu/meanz)
  # TODO for binary features check smd <- abs((p1-p2)/sqrt(p1*(1-p1)+p2*(1-p2)))
  imbalanced <- propM>maxProp | (1/maxProp)>propM  # TODO - change it to a stat test?
  if (sum(imbalanced) > 0 & verbose) {
    reportDf <- data.frame(matrix(ncol=3, nrow=sum(imbalanced)))
    rownames(reportDf) <- names(propM[imbalanced])
    colnames(reportDf) <- c('Int', 'Ext', 'Prop')
    reportDf[[1]] <- meanz[imbalanced]
    reportDf[[2]] <- mu[imbalanced]
    reportDf[[3]] <- propM[imbalanced]
    cat('Found imbalanced features:\n')
    print(reportDf)
  }
  return(imbalanced)
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
  return(list(AUC=pAuc, LogLike=pLogLike, Brier=pBrier, observedR=pMeanObservedRisk, predictedR=pMeanPredictionRisk))
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
  dbRes[['n y']] <- sum(y)
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
