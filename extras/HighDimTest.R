rm(list=ls())

library(glue)
library(LearningWithExternalStats)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./efesSimulationTests.R')
source('./plotEfesSimulationResults.R')

# Set output dir
outputDir = 'D:/projects/robustness/high-dim'  # getwd()

testParams <- getDefaultEfesTestParams(outputDir = outputDir)
bigTest <- T
bigTest2 <- F
testParams$sigma_B_Y_XA_factor <- 3
for (sigma_B_X_AH in c(0, 0.5)) {  # Degree of proximity assumption violation
  testParams$sigma_B_X_AH <- sigma_B_X_AH
  if (bigTest) {
    testParams$loadCached = F  # load or train from scratch
    testParams$n <- 300000
    testParams$p <- 500
    testParams$ntop <- 500  # maximum number of features to re-weight on
    testParams$outcomeOffset <- -log(250)
    testParams$nTest <- 10
    testParams$estimationParams[[1]]$nMaxReweight <- 5000
    testParams$estimationParams[[2]]$nMaxReweight <- 5000
  } else {
    testParams$loadCached = F  # load or train from scratch
    testParams$n <- 5000
    testParams$p <- 100
    testParams$ntop <- 100
    testParams$outcomeOffset <- -log(2)
    testParams$nTest <- 5
    testParams$estimationParams[[1]]$nMaxReweight <- 1000
    testParams$estimationParams[[2]]$nMaxReweight <- 1000
    testParams$estimationParams[[2]]$nRepetitions <- 5  # number of estimation repetitions
  }
  print(testParams)
  res <- repeatedTests(testParams)
  plotHighDimResults(testParams)

  if (bigTest && bigTest2) {
    testParams$estimationParams[[1]]$nRepetitions <- 10  # number of estimation repetitions
    testParams$estimationParams[[2]]$nMaxReweight <- 10000
    print(testParams)
    res <- repeatedTests(testParams)
    plotHighDimResults(testParams)
  }
}
