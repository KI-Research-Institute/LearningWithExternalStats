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

  alphas <- c(0.01, 0.03, 0.1, 0.3, 0.5)

  nrep <- 100
  nMaxReweight <- 20000 # Need to check if more samples lead to faster convergence

  cvxAlg  <- createExternalEstimatorSettings(
    reweightAlgorithm = cvxWeightOptimizer(),
    shortName = 'CVX 5k',
    stratifiedSampling = T,
    nMaxReweight = nMaxReweight, # maximum samples in a single round of repeated estimations
    nRepetitions = nrep,
    outputDir = outputDir,
    maxCores = 15
  )

  estimationParams <- vector(mode = 'list', length = 2)
  # estimationParams[[1]] <- cvxAlg
  threshold <- 1e-4
  esti <- cvxAlg
  esti$shortName <- 'Baseline'
  esti$reweightAlgorithm <-
    seTunedWeightOptimizer(
      outcomeCol = 'Y', alphas = alphas, outputDir=outputDir, improveTh = threshold, maxErr = 1, nIter = 500,
      nTuneIter=25)
  estimationParams[[1]] <- esti
  esti <- cvxAlg
  esti$shortName <- 'Warm 10/100'
  esti$reweightAlgorithm <- seTunedWeightOptimizer(
    outcomeCol = 'Y', alphas = alphas, outputDir=outputDir, improveTh = threshold, maxErr = 1, nIter = 500,
    nTuneIter=25)
  esti$warmStartAlgorithm <- seTunedWeightOptimizer(
    outcomeCol = 'Y', alphas = alphas, outputDir=outputDir, improveTh = threshold, maxErr = 1, nIter = 100,
    nTuneIter=5)
  estimationParams[[2]] <- esti

  minW <- 0  # minimum weight
  maxProp = 100  # TODO refine this test
  testParams <- list(
    n = NA,  # number of samples in train and tests sets
    p = NA,  # number of features
    binary = T, # type of covariates
    ntop = 1000,  # number of features used in estimation of external performance
    nTest = nTest,  # number of tests
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

proximityAndOffsetTests(cfg, testParams, nTest, loadCached)  # TODO why is there a duplication in nTest?
