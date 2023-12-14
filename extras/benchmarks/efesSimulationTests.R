#' Estimation from external statistics (EFES) simulation tests
#'


library(LearningWithExternalStats)
library(glue)
library(glmnet)
library(pROC)

source('../mlWrappers/wglmnet.R')
source('../mlWrappers/wxgboost.R')
source('../mlWrappers/wkeras.R')
source('../mlWrappers/wrappedml.R')
source('../simulations/anchorModelSimulator.R')


testSimulatedData <- function(testParams, testNum) {
  # Generate a single test set for different configuration
  # TODO change the parameters structure
  dataFileName <- file.path(
    testParams$dataDir, glue('data_{getTestName(testParams)}_{testNum}.rds'))
  cat('data file name:', dataFileName, '\n')

  if (testParams$loadCached) {
    data <- readRDS(file=dataFileName)
    trainingTime <- data$trainingTime
    d <- data$d
    model1 <- data$model1
    xFeatures <- data$xFeatures
  }
  else {
    d <- generateSimulatedData(testParams)
    # Train a model
    xFeatures <- colnames(d$internalTest)[1:(ncol(d$internalTest)-1)]
    cat('Generated simulated data, with the following means:\n')
    print(sort(colMeans(d$internalTest))[1:min(testParams$p, 10)])
    trainingStartTime <- Sys.time()
    #model1 <- cv.glmnet(sapply(d$internalTrain[xFeatures], as.numeric), d$internalTrain[['Y']],
    #                    family = "binomial", type.measure = "auc", alpha = 0.5)
    cat('Training using ', testParams$trainer$name, '\n')
    model1 <- wfit(testParams$trainer, sapply(d$internalTrain[xFeatures], as.numeric), d$internalTrain[['Y']])
    trainingTime <- Sys.time() - trainingStartTime
    saveRDS(list(d=d, model1=model1, trainingTime=trainingTime, xFeatures=xFeatures), dataFileName)
  }
  # vars1 <- wimportantGlmnet(model1, testParams$ntop)
  vars1 <- wimportant(model1)
  cat('Number of selected variables' ,length(vars1), '\n')
  wprintGlmnet(model1)

  # Predict the label probabilities in the internal test set
  internalX <- sapply(d$internalTest[xFeatures], as.numeric)
  # pInternal <- predict(model1, internalX, type = "response", s = "lambda.1se")[,1]
  pInternal <- wpredict(model1, internalX)
  internalAUC <- auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T))
  internalBrier <- mean((d$internalTest[['Y']]-pInternal)^2)
  internalCalibrationObserved <- mean(d$internalTest[['Y']])
  internalCalibrationPredicted <- mean(pInternal)

  # External results
  xExternal <- sapply(d$externalTest[xFeatures], as.numeric)
  # pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
  pExternal <- wpredict(model1, xExternal)
  extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))
  extBrier <- mean((d$externalTest[['Y']]-pExternal)^2)
  extCalibrationObserved <- mean(d$externalTest[['Y']])
  extCalibrationPredicted <- mean(pExternal)

  internalY <- d$internalTest[['Y']]
  results <- list(
    'n' = length(internalY),
    'n outcome' = sum(internalY),
    'Internal AUC' = internalAUC,
    'Internal Brier' = internalBrier,
    'Internal Calibration prediction' = internalCalibrationPredicted,
    'Internal Calibration observed' = internalCalibrationObserved,
    'External AUC' = extAuc,
    'External Brier' = extBrier,
    'External Calibration prediction' = extCalibrationPredicted,
    'External Calibration observed' = extCalibrationObserved,
    'Training Time' = trainingTime,
    'n Reweight Vars' = length(vars1)
  )

  abbrevations <- list(
    'AUC' = 'AUROC',
    'Brier' = 'Brier score',
    'Calibration prediction' = 'Global calibration mean prediction',
    'Calibration observed' = 'Global calibration observed risk'
  )
  for (i in 1:length(testParams$estimationParams)) {

    estResults <- estimatePerformance(testParams, i, d, pInternal, vars1)
    fieldNames <- sapply(getPreDiagnosticsFieldNames(), function(x) glue('{x} {i}'))
    results[fieldNames] <- estResults$results[getPreDiagnosticsFieldNames(), 'value']

    if (!is.null(estResults$weightingResults)) {
      for (metric in c('Opt err', 'n iter')) {
        a <- metric
        results[glue('Est. {metric} {i}')] <- estResults$results[a, 'value']
        results[glue('Est. {metric} {i} low')] <- estResults$results[glue('95% lower {a}'), 'value']
        results[glue('Est. {metric} {i} high')] <- estResults$results[glue('95% upper {a}'), 'value']
      }
      results[glue('Estimation Time {i}')] <- estResults$estimationTime
    }
    else {
      for (metric in c('Opt err', 'n iter')) {
        results[glue('Est. {metric} {i}')] <- NA
        results[glue('Est. {metric} {i} low')] <- NA
        results[glue('Est. {metric} {i} high')] <- NA
      }
      results[glue('Estimation Time {i}')] <- estResults$estimationTime
      cat('Failed test', i, '\n')
    }

    if (!is.null(estResults$estimation)) {
      for (metric in names(abbrevations)) {
        a <- abbrevations[[metric]]
        results[glue('Est. {metric} {i}')] <- estResults$results[a, 'value']
        results[glue('Est. {metric} {i} low')] <- estResults$results[glue('95% lower {a}'), 'value']
        results[glue('Est. {metric} {i} high')] <- estResults$results[glue('95% upper {a}'), 'value']
      }
      for (metric in c('Max Weighted SMD'))  # , 'chi2 to uniform', 'kl'
        results[glue('{metric} {i}')] <- estResults$results[metric, 'value']
    }
    else {
      for (metric in names(abbrevations)) {
        results[glue('Est. {metric} {i}')] <- NA
        results[glue('Est. {metric} {i} low')] <- NA
        results[glue('Est. {metric} {i} high')] <- NA
      }
      for (metric in c('Max Weighted SMD')) # , 'chi2 to uniform', 'kl'
        results[glue('{metric} {i}')] <- NA
      cat('Failed test', i, '\n')
    }
  }

  return(results)
}


estimatePerformance <- function(testParams, i, d, pInternal, vars1) {
  ## Estimation using reweighing
  # TODO HIGHDIM
  estimationParams <- testParams$estimationParams[[i]]

  internalK <- d$internalTest[, c(vars1, 'Y')]
  externalK <- d$externalTest[, c(vars1, 'Y')]
  transformType <- estimationParams[[1]]

  cat('\n\n--- estimating performance ---\n\n')

  dTransformedInt <- transformClassifierData(internalK, transformType, interactionVar = vars1[1])
  dTransformedExt <- transformClassifierData(externalK, transformType, interactionVar = vars1[1])

  print(colnames(dTransformedInt))
  muExt <- colMeans(dTransformedExt)
  internalData <- list(z=dTransformedInt, p = pInternal, y = internalK[['Y']])

  saveDataCSV <- F
  if (saveDataCSV) {
    csvPrefix <- glue('n{nrow(internalData$z)} m{ncol(internalData$z)} vi{testParams$sigma_B_X_AH}')
    write.csv(
      internalData$z, file = file.path(testParams$outputDir, glue('{csvPrefix} internal features.csv')))
    write.csv(
      data.frame(mean.value=muExt), file = file.path(testParams$outputDir, glue('{csvPrefix} external means.csv')))
  }

  estimatorSettings  <- createExternalEstimatorSettings(
    reweightAlgorithm = estimationParams[[2]],
    nRepetitions = testParams$nRepetitions,
    outputDir = outputDir,
    maxCores = testParams$maxCores)

  res <- estimateExternalPerformanceFromStatistics(
    internalData = internalData,
    externalStats = muExt,
    externalEstimatorSettings = estimatorSettings)
  return(res)
}


repeatedTests <- function(params) {
  testName <- getTestName(params)
  cat("Running", testName, "\n")
  r <- testSimulatedData(params, 1)
  res <- data.frame(matrix(ncol=length(r), nrow=params$nTest))
  colnames(res) <- names(r)
  res[1, ] <- r
  write.csv(res[1, ], file.path(params$outputDir, glue('{testName} 1-1.csv')))
  for (i in 2:params$nTest) {
    r <- testSimulatedData(params, i)
    res[i, ] <- r
    write.csv(res[1:i, ], file.path(params$outputDir, glue('{testName} 1-{i}.csv')))
    cat(glue('Completed test {i}'), '\n')
  }
  res['diff'] <- abs(res['Internal AUC'] - res['External AUC'])
  for (k in 1:length(params$estimationParams)) {
    res[glue('err {k}')] <- abs(sapply(res[glue('Est. AUC {k}')], as.numeric) - res['External AUC'])
  }
  write.csv(res, file.path(params$outputDir, glue('{testName}.csv')))
  return(res)
}


# getDataName <- function(params) {
#   dataName <- glue("n{params$n}-p{params$p}-top{params$ntop}-o{exp(params$outcomeOffset)}")
#   dataName <- glue("{dataName}-vi{params$sigma_B_X_AH}")
#   return(dataName)
#}


getTestName <- function(params) {
  testName <- getModelName(params)
  testName <- glue("{testName}-{params$trainer$name}-vi{params$sigma_B_X_AH}")
  return(testName)
}
