library(glmnet)


predictGLM <- function(xy, glModel) {
  xFeatures <- colnames(xy[1:(ncol(xy)-1)])
  X <- sapply(xy[xFeatures], as.numeric)
  p <- predict(glModel, X, type = "response", s = "lambda.1se")[,1]
  return(p)
}


estimatePerformanceWithInternalData <- function(z, pInternal, muExt, wOptimizer) {

  externalEstimatorSettings <- createExternalEstimatorSettings(
    reweightAlgorithm = wOptimizer,
    nMaxReweight = nrow(z),
    nRepetitions = 2,
    maxCores = 1,
    maxWSMD = 0.2
  )
  internalData <- list(z=z, p = pInternal, y = z[['Y']])  # TODO consider changing the api

  estimatedLRResults <- estimateExternalPerformanceFromStatistics(
    internalData = internalData,
    externalStats = muExt,
    externalEstimatorSettings = externalEstimatorSettings,
    createEstimationLogger = F
  )

  return(estimatedLRResults)

}

# TODO Consider removing this
estimateExternalAUCWithInternalData <- function(z, pInternal, muExt, wOptimizer) {

  estimatedLRResults <- estimatePerformanceWithInternalData(z, pInternal, muExt, wOptimizer)
  estimatedAUC <- estimatedLRResults$estimation['AUROC', 'value']
  return(estimatedAUC)

}


replaceSimpleLoggerThreshold <- function(threshold) {
  appender <- ParallelLogger::createConsoleAppender(layout = layoutTimestamp)
  logger <- ParallelLogger::createLogger(name = "SIMPLE", threshold = threshold, appenders = list(appender))
  unregisterLogger("SIMPLE")
  registerLogger(logger = logger)
}
