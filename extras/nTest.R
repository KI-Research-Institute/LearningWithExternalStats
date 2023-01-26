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
  maxAlpha <- 1e-0
  nAlphas <- 4

  # alphas <- minAlpha*exp(log(maxAlpha/minAlpha)*(0:nAlphas)/nAlphas)
  alphas <- c(0.1, 0.2, 0.5)

  nrep <- 15

  cvxAlg  <- createExternalEstimatorSettings(
    reweightAlgorithm = cvxWeightOptimizer(),
    shortName = 'CVX 5k',
    stratifiedSampling = T,
    nMaxReweight = NA, # maximum samples in a single round of repeated estimations
    nRepetitions = nrep,
    outputDir = outputDir,
    maxCores = 15
  )

  threshold <- 1e-4
  estimationParams <- vector(mode = 'list', length = 2)
  esti <- cvxAlg
  esti$shortName <- '10k'
  esti$reweightAlgorithm <- seTunedWeightOptimizer(
    alphas = alphas, outputDir=outputDir, improveTh = threshold, maxErr = 1, nIter = 1000, nTuneIter=50)
  esti$nMaxReweight <- 10000
  estimationParams[[1]] <- esti
  esti <- cvxAlg
  esti$shortName <- '20k'
  esti$reweightAlgorithm <- seTunedWeightOptimizer(
    alphas = alphas, outputDir=outputDir, improveTh = threshold, maxErr = 1, nIter = 1000, nTuneIter=50)
  esti$nMaxReweight <- 20000
  estimationParams[[2]] <- esti


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
