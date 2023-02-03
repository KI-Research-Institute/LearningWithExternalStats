rm(list=ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(glue)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./offset-test.R')

# cfg <- list(p=100, n=2000, outcomeOffsets=-log(c(2, 1)))
cfg <-list(p=500, n=3e5, outcomeOffsets=-log(c(250)))
nTest <- 10
loadCached = T

getGeneralizationTestsParams <- function(outputDir) {

  nrep <- 15

  esti1  <- createExternalEstimatorSettings(
    reweightAlgorithm = seTunedWeightOptimizer(outputDir=outputDir, nIter = 5000),
    nRepetitions = nrep,
    outputDir = outputDir,
    maxCores = 15
  )


  estimationParams <- vector(mode = 'list', length = 2)
  estimationParams[[1]] <- esti1
  esti2 <- esti1
  esti2$shortName <- 'Th 5e-2'
  esti2$reweightAlgorithm <- seTunedWeightOptimizer(outputDir=outputDir, nIter = 5000, improveTh=5e-5)
  estimationParams[[2]] <- esti2

  testParams <- list(
    n = NA,  # number of samples in train and tests sets
    p = NA,  # number of features
    binary = T, # type of covariates
    ntop = 1000,  # number of features used in estimation of external performance
    nTest = nTest,  # number of tests
    # Simulation model
    outcomeOffset = NA,  # offset of the outcome logistic model, determines outcome prevalence
    sigma_B_X_AH = NA,  # degree of porximity assumption violation
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
