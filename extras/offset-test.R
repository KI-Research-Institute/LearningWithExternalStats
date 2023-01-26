rm(list=ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(glue)

# TODO handle zero weights

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./efesSimulationTests.R')
source('./plotEfesSimulationResults.R')


proximityAndOffsetTests <- function(cfg, testParams, nTest, loadCached) {
  p <- cfg$p
  n <- cfg$n #
  outcomeOffsets <- cfg$outcomeOffsets
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

    legendLabels <- vector(mode = 'list', length = k)
    for (i in 1:k)
      legendLabels[[i]] <- testParams$estimationParams[[i]]$shortName

    summaryName <- glue('n{n}p{p}s{sigma_B_X_AH}o{exp(-outcomeOffsets[1])}-{exp(-outcomeOffsets[nOffsets])}')
    write.csv(aucErrors, file.path(testParams$outputDir, glue('{summaryName} AUC.csv')))
    write.csv(brierErrors, file.path(testParams$outputDir, glue('{summaryName} Brier.csv')))

    plotOffsetResults(testParams$outputDir, summaryName, legendLabels)

  }
}
