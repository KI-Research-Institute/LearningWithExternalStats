rm(list=ls())

library(glue)
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./plotEfesSimulationResults.R')


outputDir <- 'D:/projects/robustness/offset-test-3e+05-500'

resName <- 'n3e+05-p500-top1000-o0.004-vi0.5-glmnet 1-6'
#outputDir <- 'D:/projects/robustness/offset-test-20000-1000'
#resName <- 'n2000-p100-top1000-o0.5-vi0-glmnet 1-2'
#outputDir <- 'D:/projects/robustness/offset-test-2000-100'
# resName <- 'n10000-p200-top1000-o0.1-vi0-glmnet 1-3'

legendLabels <- seq(1,1,1)
plotResults(outputDir, resName, legendLabels = legendLabels)
plotPropoertiesBox(outputDir, resName, legendLabels = legendLabels)
