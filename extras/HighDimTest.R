rm(list=ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(glue)
library(LearningWithExternalStats)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./efesSimulationTests.R')
source('./plotEfesSimulationResults.R')

# Set output dir

# In primal dual, batches should be large

p <- 20  # 10 500
n <- 3000  # 1e4 300000
outcomeOffset <- -log(4)  # 4 250
for (sigma_B_X_AH in c(0)) {  # Degree of proximity assumption violation

  testParams <- getDefaultEfesTestParams(outputDir = 'D:/projects/robustness/high-dim-small')
  testParams$sigma_B_Y_XA_factor <- 3
  testParams$sigma_B_X_AH <- sigma_B_X_AH
  testParams$loadCached = F  # load or train from scratch
  testParams$n <- n
  testParams$p <- p
  testParams$ntop <- 500  # maximum number of features to re-weight on
  testParams$outcomeOffset <- outcomeOffset
  testParams$nTest <- 3
  # testParams$estimationParams[[1]]$nTuneIter <- 10 # 10 is Too little with batch size of 1000
  testParams$estimationParams[[1]]$batchSize <- 1000000
  testParams$estimationParams[[1]]$polyBeta <- 0  # 0 is fixed step size for deterministic
  # testParams$estimationParams[[1]]$nIter <- 1000
  # testParams$estimationParams[[1]]$maxNuNorm <- 10
  # testParams$estimationParams[[1]]$improveTh <- -1
  # testParams$estimationParams[[1]]$nProbe <- 50

  res <- repeatedTests(testParams)
  plotHighDimResults(testParams)
}
