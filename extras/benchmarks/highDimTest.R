rm(list=ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(glue)
library(LearningWithExternalStats)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./offset-test.R')
source('../optimize/kerasWeightOptimizer.R')

testLevel <- 'a'

# Parameters

loadCached <- F

testParams <- list(
  binary = T, # type of covariates
  p=NaN,  # 2
  n=NaN,  # 2000
  outcomeOffset = NaN, # -log(2),
  sigma_B_Y_X_factor = 3,  # Testing smaller factor means less predictive power of variables ...
  sigma_B_Y_XA_factor = 3,  # Testing ...
  envOffset = 1,  # Testing ...
  loadCached = loadCached,  # load or train from scratch
  trainer = wglmnet(), # wXGBoost()
  nRepetitions = 10,
  maxCores = 10
)

switch (testLevel,
  a = {
    testParams[c('p', 'n_binary', 'n', 'outcomeOffset')] <- c(2, 2, 2000, -log(2))
    nTest <- 2  # Should be >1
  },
  b = {
    testParams[c('p','n_binary', 'n', 'outcomeOffset')] <- c(10, 10, 10000, -log(5))
    nTest <- 10  # Should be >1
  },
  c = {
    testParams[c('p', 'n_binary', 'n','outcomeOffset')] <- c(500, 498, 3e5, -log(250))
    nTest <- 10  # Should be >1
  }
)


modelName <- getModelName(testParams)
outputDir <- file.path('../../../../../projects/robustness/bench11', modelName)
dataDir <- file.path('../../../../../projects/robustness/evaluationData', modelName)


reweightEg500 <- seTunedWeightOptimizer(outcomeCol = 'Y', outputDir=outputDir, nIter = 500)
reweightEg5000 <- seTunedWeightOptimizer(outcomeCol = 'Y', outputDir=outputDir, nIter = 5000)
# reweightKerasId <- kerasWeightOptimizer(optimizer = 'adam', nIter = 5000, parameterization = 'id')
# reweightKerasLog <- kerasWeightOptimizer(optimizer = 'adam', nIter = 5000, parameterization = 'log')

# Every item includes: transformation, reweight-algorithm and name
estimationParams <- list(
  list('Current flat', reweightEg500, 'C 500'),
  list('New flat',  reweightEg5000, 'EG N'),
  # list('New flat',  reweightKerasId, 'KerasId N'),
  # list('New flat',  reweightKerasLog, 'KerasLog N'),
  list('Interaction v0',  reweightEg5000, 'EG I'), # ,
  list('Interaction',  reweightEg5000, 'EG I') # ,
  # list('Interaction',  reweightKerasId, 'KerasId I'),
  # list('Interaction',  reweightKerasLog, 'KerasLog I')
)

proximityAndOffsetTests(testParams, estimationParams, nTest, loadCached, outputDir, dataDir)
