#' Estimation from external statistics (EFES) simulation tests
#'


library(LearningWithExternalStats)
library(glue)
library(glmnet)
library(pROC)

source('./mlWrappers/wglmnet.R')
source('./mlWrappers/wxgboost.R')
source('./mlWrappers/wrpart.R')
# source('./mlWrappers/wkeras.R')
source('./mlWrappers/wrappedml.R')
source('./simulations/anchorModelSimulator.R')


getDefaultEfesTestParams <- function(outputDir) {

  minEpsilon0 <- 1e-4
  maxEpsilon0 <- 1e-0
  nEpsilon0 <- 4

  minAlpha <- 0.1
  maxAlpha <- 10
  nAlphas <- 4

  alphas <- minAlpha*exp(log(maxAlpha/minAlpha)*(0:nAlphas)/nAlphas)
  epsilon0s = minEpsilon0*exp(log(maxEpsilon0/minEpsilon0)*(0:nEpsilon0)/nEpsilon0)
  epsilon0sDeterministic = minEpsilon0*exp(log(1e-1/minEpsilon0)*(0:nEpsilon0)/nEpsilon0)

  # estDual <- sgdTunedWeightOptimizer(
  #   outputDir=outputDir, epsilon0s=epsilon0sDeterministic, batchSize = NA, polyBeta = 0, improveTh=1e-4, nProbe = 10)
  # estDualPrimal <- sgdTunedWeightOptimizer(
  #   outputDir=outputDir, epsilon0s=epsilon0sDeterministic, batchSize = NA, polyBeta = 0, improveTh=1e-4, nProbe = 10,
  #   primal = T)

  maxProp = 100  # TODO refine this test
  cvxAlg  <- createExternalEstimatorSettings(
    reweightAlgorithm = cvxWeightOptimizer(),
    shortName = 'CVX 5k',
    stratifiedSampling = T,
    nMaxReweight = 5000, # maximum samples in a single round of repeated estimations
    nRepetitions = 10,
    maxProp = maxProp,
    outputDir = outputDir,
    maxCores = 15
  )


  nrep <- 15
  nMaxReweight <- 500 # 5000
  # 1.
  est1 <- cvxAlg
  est1$nRepetitions = nrep
  est1$shortName <- 'W-MSE few iters'
  est1$nMaxReweight <- nMaxReweight
  est1$reweightAlgorithm <-
    seTunedWeightOptimizer(
      alphas = alphas, outputDir=outputDir, improveTh = 1e-4, maxErr = 1, nIter = 200, nTuneIter=20)

  # 2.
  est2 <- cvxAlg
  est2$nRepetitions = nrep
  est2$shortName <- 'W-MSE few warm start'
  est2$nMaxReweight <- nMaxReweight
  est2$warmStartAlgorithm <-
    seTunedWeightOptimizer(
      alphas = alphas, outputDir=outputDir, improveTh = 1e-4, maxErr = 1, nIter = 40, nTuneIter=4)
  est2$reweightAlgorithm <-
    seTunedWeightOptimizer(
      alphas = alphas, outputDir=outputDir, improveTh = 1e-4, maxErr = 1, nIter = 200, nTuneIter=20)
  # 3.
  est3 <- cvxAlg
  est3$nRepetitions = nrep
  est3$shortName <- 'W-MSE many iters'
  est3$nMaxReweight <- nMaxReweight
  est3$reweightAlgorithm <-
    seTunedWeightOptimizer(alphas = alphas, outputDir=outputDir, improveTh = 1e-4, maxErr = 1, nIter = 500, nTuneIter=50)

  # 4.
  est4 <- cvxAlg
  est4$nRepetitions = nrep
  est4$shortName <- 'W-MSE many warm start'
  est4$nMaxReweight <- nMaxReweight
  est4$warmStartAlgorithm <-
    seTunedWeightOptimizer(
      alphas = alphas, outputDir=outputDir, improveTh = 1e-4, maxErr = 1, nIter = 40, nTuneIter=4)
    # sgdTunedWeightOptimizer(outputDir=outputDir, epsilon0s=epsilon0s, batchSize = 1000, polyBeta = 1)
  est4$reweightAlgorithm <-
    seTunedWeightOptimizer(
      alphas = alphas, outputDir=outputDir, improveTh = 1e-4, maxErr = 1, nIter = 500, nTuneIter=50)


  minW <- 0  # minimum weight
  maxProp = 100  # TODO refine this test
  testParams <- list(
    n = 3e5,  # number of samples in train and tests sets
    p = 500,  # number of features
    binary = T, # type of covariates
    ntop = 1000,  # number of features used in estimation of external performance
    nTest = 20,  # number of tests
    # Simulation model
    outcomeOffset = -log(250),  # offset of the outcome logistic model, determines outcome prevalence
    sigma_B_X_AH = 0,  # degree of porximity assumption violation
    sigma_B_Y_X_factor = 4,
    sigma_B_Y_XA_factor = 4,
    loadCached = F,  # load or train from scratch
    envOffset = 5,
    # Estimation model
    trainer = wglmnet(),  # wXGBoost()
    outputDir = outputDir,
    # Reweighing parameters
    estimationParams = list(
      #cvxAlg,
      est1,
      est2,
      est3,
      est4
    )
  )
  # Notes: SCS failed with n=1000 p=100 in 1/2 experiments. n=5000, p=10 5/7
  # ECOS_BB quite similar to ECOS (BB is branch and bound and suitable for mixed integer problems)
  return(testParams)
}






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
    print(sort(colMeans(d$internalTest))[floor((0:10)/10*(testParams$p-1)+1)])
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
    'Calibration observed' = 'Global calibration observed risk',
    'Opt err' = 'Opt err',
    'n iter' = 'n iter'
  )
  for (i in 1:length(testParams$estimationParams)) {
    estResults <- estimatePerformance(testParams, i, d, pInternal, vars1)

    fieldNames <- sapply(getPreDiagnosticsFieldNames(), function(x) glue('{x} {i}'))
    results[fieldNames] <- estResults$results[getPreDiagnosticsFieldNames(), 'value']

    # if (estResults$status == 'Success') {
    if (!is.null(estResults$estimation)) {
      for (metric in names(abbrevations)) {
        a <- abbrevations[[metric]]
        results[glue('Est. {metric} {i}')] <- estResults$results[a, 'value']
        results[glue('Est. {metric} {i} low')] <- estResults$results[glue('95% lower {a}'), 'value']
        results[glue('Est. {metric} {i} high')] <- estResults$results[glue('95% upper {a}'), 'value']
      }
      results[glue('Estimation Time {i}')] <- estResults$estimationTime
      for (metric in c('Max Weighted SMD'))  # , 'chi2 to uniform', 'kl'
        results[glue('{metric} {i}')] <- estResults$results[metric, 'value']
    }
    else {
      for (metric in names(abbrevations)) {
        results[glue('Est. {metric} {i}')] <- NA
        results[glue('Est. {metric} {i} low')] <- NA
        results[glue('Est. {metric} {i} high')] <- NA
      }
      results[glue('Estimation Time {i}')] <- estResults$estimationTime
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

  internalK <- d$internalTest[, c(vars1, 'Y')]
  externalK <- d$externalTest[, c(vars1, 'Y')]
  dTransformedInt <- computeTable1LikeTransformation(internalK, outcomeBalance=TRUE)
  dTransformedExt <- computeTable1LikeTransformation(externalK, outcomeBalance=TRUE)
  muExt <- colMeans(dTransformedExt)
  internalData <- list(z=dTransformedInt, p = pInternal, y = internalK[['Y']])
  estimationParams <- testParams$estimationParams[[i]]

  saveDataCSV <- F
  if (saveDataCSV) {
    csvPrefix <- glue('n{nrow(internalData$z)} m{ncol(internalData$z)} vi{testParams$sigma_B_X_AH}')
    write.csv(
      internalData$z, file = file.path(testParams$outputDir, glue('{csvPrefix} internal features.csv')))
    write.csv(
      data.frame(mean.value=muExt), file = file.path(testParams$outputDir, glue('{csvPrefix} external means.csv')))
  }

  res <- estimateExternalPerformanceFromStatistics(
    internalData = internalData,
    externalStats = muExt,
    externalEstimatorSettings = estimationParams)
  return(res)
}


repeatedTests <- function(params) {
  testName <- getTestName(params)
  cat("Running", testName, "\n")
  r <- testSimulatedData(params, 1)
  res <- data.frame(matrix(ncol=length(r), nrow=params$nTest))
  colnames(res) <- names(r)
  res[1, ] <- r
  print(res)
  print(file.path(params$outputDir, glue('{testName} 1-1.csv')))
  write.csv(res[1, ], file.path(params$outputDir, glue('{testName} 1-1.csv')))
  for (i in 2:params$nTest) {
    r <- testSimulatedData(params, i)
    res[i, ] <- r
    # print(res)
    write.csv(res[1:i, ], file.path(params$outputDir, glue('{testName} 1-{i}.csv')))
    cat(glue('Completed test {i}'), '\n')
  }
  res['diff'] <- abs(res['Internal AUC'] - res['External AUC'])
  for (k in 1:length(params$estimationParams)) {
    res[glue('err {k}')] <- abs(sapply(res[glue('Est. AUC {k}')], as.numeric) - res['External AUC'])
  }
  print(res)
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
