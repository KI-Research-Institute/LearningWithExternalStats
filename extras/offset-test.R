rm(list=ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(glue)
library(LearningWithExternalStats)
# library(abind)

# TODO handle zero weights

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./efesSimulationTests.R')
source('./plotEfesSimulationResults.R')


# Set output dir

# In primal dual, batches should be large

configs <- list(
  list(p=100, n=2000, outcomeOffsets=-log(c(2, 1))),  # TODO debug this
  list(p=200, n=10000, outcomeOffsets=-log(c(10, 3))),
  list(p=1000, n=20000, outcomeOffsets=-log(c(20, 5))),
  list(p=500, n=3e5, outcomeOffsets=-log(c(250, 16)))
)

cfg <- configs[[1]]
p <- cfg$p
n <- cfg$n #
outcomeOffsets <- cfg$outcomeOffsets

nTest <- 2
loadCached = T

testParams <- getDefaultEfesTestParams(outputDir = glue('D:/projects/robustness/offset-test-{n}-{p}'))
nOffsets <- length(outcomeOffsets)

dir.create(testParams$outputDir, showWarnings = F)


for (sigma_B_X_AH in c(0, 0.5)) {  # Degree of proximity assumption violation

  aucErrors  <- array(dim=c(nOffsets * length(testParams$estimationParams), nTest))
  brierErrors  <- array(dim=c(nOffsets * length(testParams$estimationParams), nTest))

  for (i in 1:nOffsets) {
    outcomeOffset <- outcomeOffsets[i]
    testParams$sigma_B_Y_XA_factor <- 4
    testParams$sigma_B_X_AH <- sigma_B_X_AH
    testParams$loadCached = loadCached  # load or train from scratch
    testParams$n <- n
    testParams$p <- p
    testParams$outcomeOffset <- outcomeOffset
    testParams$nTest <- nTest

    res <- repeatedTests(testParams)
    plotHighDimResults(testParams)

    k <- length(testParams$estimationParams)
    for (j in 1:k) {
      aucErrors[(i-1)* length(testParams$estimationParams)+j, ] <-
        res[, glue('Est. AUC {j}')] - res[, 'External AUC']  # error1
      brierErrors[(i-1)* length(testParams$estimationParams)+j, ] <-
        res[, glue('Est. Brier {j}')] - res[, 'External Brier']  # error1
    }

  }

  colors <- c('red', 'purple', 'orange', 'blue', 'brown')
  legendLabels <-  c('Dual', 'SE sum', 'CVX')  #  ,

  legendLabels <- vector(mode = 'list', length = k)
  for (i in 1:k)
    legendLabels[[i]] <- testParams$estimationParams[[i]]$shortName


  logName <- glue('n{n}p{p}s{sigma_B_X_AH}o{exp(-outcomeOffset)} AUC.png')
  png(filename = file.path(testParams$outputDir, logName), width = 480, height = 480)
  boxplot(t(aucErrors), col = colors[1:k])
  legend('bottomright', fill = colors[1:k], legend=legendLabels, cex=1, inset=0.02)
  dev.off()
  logName <- glue('n{n}p{p}s{sigma_B_X_AH}o{exp(-outcomeOffset)} Brier.png')
  png(filename = file.path(testParams$outputDir, logName), width = 480, height = 480)
  boxplot(t(brierErrors), col = colors[1:k])
  legend('bottomright', fill = colors[1:k], legend=legendLabels, cex=1, inset=0.02)
  dev.off()

}
