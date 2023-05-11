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
nTest <- 10  # Should be >1
loadCached = T

testParams <- list(
  binary = T, # type of covariates
  p=NaN,  # 2
  n=NaN,  # 2000
  outcomeOffset = NaN, # -log(2),
  sigma_B_Y_X_factor = 3,  # Testing smaller factor means less predictive power of variables ...
  sigma_B_Y_XA_factor = 3,  # Testing ...
  envOffset = 1,  # Testing ...
  loadCached = loadCached,  # load or train from scratch
  trainer = wglmnet()  # wXGBoost()
)
# testParams[c('p', 'n', 'outcomeOffset')] <- c(2, 2000, -log(2))
testParams[c('p', 'n', 'outcomeOffset')] <- c(500, 3e5, -log(250))

modelName <- getModelName(testParams)
outputDir <- file.path('D:/projects/robustness/bench02', modelName)
dataDir <- file.path('D:/projects/robustness/evaluationData', modelName)

esti  <- createExternalEstimatorSettings(
  reweightAlgorithm = seTunedWeightOptimizer(outputDir=outputDir),
  nRepetitions = 15,
  outputDir = outputDir,
  maxCores = 15
)
estimationParams <- vector(mode = 'list', length = 1)
estimationParams[[1]] <- esti

proximityAndOffsetTests(testParams, estimationParams, nTest, loadCached, outputDir, dataDir)
