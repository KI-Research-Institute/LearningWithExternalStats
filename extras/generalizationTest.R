rm(list=ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(glue)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./offset-test.R')

# cfg <- list(p=100, n=2000, outcomeOffsets=-log(c(2, 1)))
cfg <-list(p=500, n=3e5, outcomeOffsets=-log(c(250, 16)))
nTest <- 10
loadCached = T

getGeneralizationTestsParams <- function(outputDir) {

  minAlpha <- 1e-2
  maxAlpha <- 1e-1
  nAlphas <- 4

  alphas <- minAlpha*exp(log(maxAlpha/minAlpha)*(0:nAlphas)/nAlphas)

  nrep <- 30
  nMaxReweight <- 20000 # Need to check if more samples lead to faster convergence

  maxProp = 100  # TODO refine this test
  cvxAlg  <- createExternalEstimatorSettings(
    reweightAlgorithm = cvxWeightOptimizer(),
    shortName = 'CVX 5k',
    stratifiedSampling = T,
    nMaxReweight = nMaxReweight, # maximum samples in a single round of repeated estimations
    nRepetitions = nrep,
    maxProp = maxProp,
    outputDir = outputDir,
    maxCores = 15
  )

  thresholds <- 10**seq(-4, -4,-0.5)
  estimationParams <- vector(mode = 'list', length = 3*length(thresholds))
  for (i in 1:length(thresholds)) {
    esti <- cvxAlg
    esti$shortName <- glue('Th{log10(thresholds[i])}m{esti$reweightAlgorithm$momentumC}')
    esti$reweightAlgorithm <-
      seTunedWeightOptimizer(
        outcomeCol = 'Y', alphas = alphas, outputDir=outputDir, improveTh = thresholds[i], maxErr = 1, nIter = 500,
        nTuneIter=50, momentumC = 0.9)
    estimationParams[[3*i-2]] <- esti
    esti$reweightAlgorithm$momentumC <- 0.8
    esti$shortName <- glue('Th{log10(thresholds[i])}m{esti$reweightAlgorithm$momentumC}')
    estimationParams[[3*i-1]] <- esti
    esti$reweightAlgorithm$momentumC <- 0.5
    # esti$nMaxReweight <- 50000
    esti$shortName <- glue('Th{log10(thresholds[i])}m{esti$reweightAlgorithm$momentumC}')
    estimationParams[[3*i]] <- esti
  }

  minW <- 0  # minimum weight
  maxProp = 100  # TODO refine this test
  testParams <- list(
    n = NA,  # number of samples in train and tests sets
    p = NA,  # number of features
    binary = T, # type of covariates
    ntop = 1000,  # number of features used in estimation of external performance
    nTest = 20,  # number of tests
    # Simulation model
    outcomeOffset = NA,  # offset of the outcome logistic model, determines outcome prevalence
    sigma_B_X_AH = 0,  # degree of porximity assumption violation
    sigma_B_Y_X_factor = 4,
    sigma_B_Y_XA_factor = 4,
    loadCached = F,  # load or train from scratch
    envOffset = 5,
    # Estimation model
    trainer = wglmnet(),  # wXGBoost()
    outputDir = outputDir,
    # Reweighing parameters
    estimationParams = estimationParams
  )
  return(testParams)
}


testParams <- getGeneralizationTestsParams(outputDir = glue('D:/projects/robustness/offset-test-{cfg$n}-{cfg$p}'))

proximityAndOffsetTests(cfg, testParams, nTest, loadCached)
