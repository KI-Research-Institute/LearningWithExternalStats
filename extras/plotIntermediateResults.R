rm(list=ls())

library(glue)
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./plotEfesSimulationResults.R')

outputDir <- 'D:/projects/robustness/bench01/p500-n3e+05-bTRUE-o250-fX3-fXA3-e1'
resName <- 'p500-n3e+05-bTRUE-o250-fX3-fXA3-e1-glmnet-vi0 1-4'

legendLabels <- seq(1,3,1)
plotResults(outputDir, resName, legendLabels = legendLabels)
plotPropoertiesBox(outputDir, resName, legendLabels = legendLabels)
