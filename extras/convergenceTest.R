rm(list=ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(glue)
library(LearningWithExternalStats)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./offset-test.R')

# Parameters
nTest <- 10
loadCached = F

testParams <- list(
  binary = T, # type of covariates
  p=500,
  n=3e5,
  outcomeOffset = -log(250),
  sigma_B_Y_X_factor = 3,  # Testing smaller factor means less predictive power of variables ...
  sigma_B_Y_XA_factor = 3,  # Testing ...
  envOffset = 1,  # Testing ...
  loadCached = loadCached,  # load or train from scratch
  trainer = wglmnet()  # wXGBoost()
)
# modelParams[c('p', 'n', 'outcomeOffset')] <- c(2, 3e-5, -log(250))

modelName <- getModelName(testParams)
outputDir <- file.path('D:/projects/robustness/bench01', modelName)
dataDir <- file.path('D:/projects/robustness/evaluationData', modelName)

esti  <- createExternalEstimatorSettings(
  reweightAlgorithm = seTunedWeightOptimizer(outputDir=outputDir),
  nRepetitions = 15,
  outputDir = outputDir,
  maxCores = 15
)
estimationParams <- vector(mode = 'list', length = 3)
estimationParams[[1]] <- esti
esti$shortName <- 'previous'
esti$reweightAlgorithm$absTol <- 1e-4
estimationParams[[2]] <- esti
estimationParams[[3]] <- createExternalEstimatorSettings(
  reweightAlgorithm = cvxWeightOptimizer(),
  nRepetitions = esti$nRepetitions,
  outputDir = outputDir,
  maxCores = esti$maxCores
)

proximityAndOffsetTests(testParams, estimationParams, nTest, loadCached, outputDir, dataDir)
