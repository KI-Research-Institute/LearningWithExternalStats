rm(list=ls())

library(LearningWithExternalStats)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./efestSimulationTests.R')


testBinarySimulation <- function() {
  h <- getDefaultModelHyperParams()
  h$binary <- T
  m <- modelParams(h)
  m$A <- 0
  d <- sampleModel(m, 1000)
  xFeatures <- paste('X', 1:h$p, sep='')
  X <- d[xFeatures]
  plot(colMeans(X), m$freqs)
}


testBinarySimulation()
