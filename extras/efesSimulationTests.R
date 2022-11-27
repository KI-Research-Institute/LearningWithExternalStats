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
  nmax <- 5000  # maximum samples in a single round of repeated estimations
  minW <- 0  # minimum weight
  outputDir = outputDir
  maxProp = 100  # TODO refine this test
  testParams <- list(
    n = 3e5,  # number of samples in train and tests sets
    p = 500,  # number of features
    binary = T, # type of covariates
    ntop = 500,  # number of features used in estimation of external performance
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
      createExternalEstimatorSettings(
        divergence = 'entropy',  # entropy, chi2
        lambda = 1e-1,
        minW = minW,
        optimizationMethod = 'dual',  # dual, primal
        nMaxReweight = nmax,
        stratifiedSampling = T,
        nRepetitions = 1,
        maxProp = maxProp,
        outputDir = outputDir,
        maxCores = 15
      ),
      createExternalEstimatorSettings(
        divergence = 'entropy',  # entropy, chi2
        lambda = 1e-1,
        minW = minW,
        optimizationMethod = 'dual',  # dual, primal
        stratifiedSampling = T,
        nMaxReweight = nmax,
        nRepetitions = 5,
        maxProp = maxProp,
        outputDir = outputDir,
        maxCores = 15
      )
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
    testParams$outputDir, glue('data_{getDataName(testParams)}_{testParams$trainer$name}_{testNum}.rds'))
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

  # External results
  xExternal <- sapply(d$externalTest[xFeatures], as.numeric)
  # pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
  pExternal <- wpredict(model1, xExternal)
  extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))

  observedRisk <- mean(d$externalTest[['Y']])
  predictedRisk <- mean(pExternal)
  cat('------- Observerved Risk', observedRisk, '----------\n')
  cat('------- Predicted   Risk', predictedRisk, '----------\n')

  internalY <- d$internalTest[['Y']]
  results <- list(
    'n' = length(internalY),
    'n outcome' = sum(internalY),
    'Internal AUC' = internalAUC,
    'External AUC' = extAuc,
    'Training Time' = trainingTime,
    'n Reweight Vars' = length(vars1)
  )

  for (i in 1:length(testParams$estimationParams)) {
    estResults <- estimatePerformance(testParams, i, d, pInternal, vars1)
    if (estResults$status == 'Success') {
      results[glue('Est. Ext. AUC {i}')] <- estResults$estimation['AUROC', 'value']
      results[glue('Estimation Time {i}')] <- estResults$estimationTime
      for (metric in c('Max Weighted SMD', 'chi2 to uniform', 'kl'))
        results[glue('{metric} {i}')] <- estResults$estimation[metric, 'value']
    }
    else {
      results[glue('Est. Ext. AUC {i}')] <- NaN
      results[glue('Estimation Time {i}')] <- NaN
      for (metric in c('Max Weighted SMD', 'chi2 to uniform', 'kl'))
        results[glue('{metric} {i}')] <- NaN
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
  write.csv(res[1, ], file.path(params$outputDir, glue('{testName} 1-1.csv')))
  for (i in 2:params$nTest) {
    r <- testSimulatedData(params, i)
    res[i, ] <- r
    print(res)
    write.csv(res[1:i, ], file.path(params$outputDir, glue('{testName} 1-{i}.csv')))
  }
  res['diff'] <- abs(res['Internal AUC'] - res['External AUC'])
  for (k in 1:length(params$estimationParams))
    res[glue('err {k}')] <- abs(res[glue('Est. Ext. AUC {k}')] - res['External AUC'])
  print(res)
  write.csv(res, file.path(params$outputDir, glue('{testName}.csv')))
  return(res)
}


getDataName <- function(params) {
  dataName <- glue("n{params$n}-p{params$p}-top{params$ntop}-o{exp(params$outcomeOffset)}")
  dataName <- glue("{dataName}-vi{params$sigma_B_X_AH}")
  return(dataName)
}


getTestName <- function(params) {
  testName <- getDataName(params)
  testName <- glue("{testName}-{params$estimationParams[[1]]$nMaxReweight}-{params$estimationParams[[1]]$nRepetitions}")
  testName <- glue("{testName}-{params$estimationParams[[2]]$nMaxReweight}-{params$estimationParams[[2]]$nRepetitions}")
  testName <- glue("{testName}-{params$trainer$name}")
  return(testName)
}
