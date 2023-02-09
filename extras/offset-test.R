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


getModelName <- function(p) {
  return(glue('p{p$p}-n{p$n}-b{p$binary}-o{exp(-p$outcomeOffset)}-fX{p$sigma_B_Y_X_factor}-fXA{p$sigma_B_Y_XA_factor}-e{p$envOffset}'))
}


proximityAndOffsetTests <- function(testParams, estimationParams, nTest, loadCached, outputDir, dataDir) {  # cfg, testParams
  p <- testParams$p
  n <- testParams$n #
  outcomeOffsets <- c(testParams$outcomeOffset) # TODO currently disabled multiple offsets
  nOffsets <- length(outcomeOffsets)
  dir.create(outputDir, showWarnings = F)
  dir.create(dataDir, showWarnings = F)

  for (sigma_B_X_AH in c(0, 0.5)) {  # Degree of proximity assumption violation

    aucErrors  <- array(dim=c(nOffsets * length(estimationParams), nTest))
    brierErrors  <- array(dim=c(nOffsets * length(estimationParams), nTest))

    for (i in 1:nOffsets) {
      outcomeOffset <- outcomeOffsets[i]
      testParams$sigma_B_X_AH <- sigma_B_X_AH
      testParams$loadCached = loadCached  # load or train from scratch
      testParams$n <- n
      testParams$p <- p
      testParams$outcomeOffset <- outcomeOffset
      testParams$nTest <- nTest
      testParams$dataDir <- dataDir
      testParams$outputDir <- outputDir
      testParams$ntop <- 1000
      testParams$estimationParams <- estimationParams  # TODO GIVE THIS TO THE FUNCTION

      res <- repeatedTests(testParams)
      plotHighDimResults(testParams)

      k <- length(estimationParams)
      for (j in 1:k) {
        aucErrors[(i-1)* length(estimationParams)+j, ] <-
          res[, glue('Est. AUC {j}')] - res[, 'External AUC']  # error1
        brierErrors[(i-1)* length(estimationParams)+j, ] <-
          res[, glue('Est. Brier {j}')] - res[, 'External Brier']  # error1
      }

    }

    legendLabels <- vector(mode = 'list', length = k)
    for (i in 1:k)
      legendLabels[[i]] <- estimationParams[[i]]$shortName

    summaryName <- glue('{getModelName(testParams)}-s{sigma_B_X_AH}')
    write.csv(aucErrors, file.path(outputDir, glue('{summaryName} AUC.csv')))
    write.csv(brierErrors, file.path(outputDir, glue('{summaryName} Brier.csv')))

    plotOffsetResults(outputDir, summaryName, legendLabels)

  }
}
