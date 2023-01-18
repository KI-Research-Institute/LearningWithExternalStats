rm(list=ls())

library(glue)
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./plotEfesSimulationResults.R')


outputDir <- 'D:/projects/robustness/offset-test-3e+05-500'
resName <- 'n3e+05-p500-top1000-o0.004-vi0-glmnet 1-3'
# outputDir <- 'D:/projects/robustness/offset-test-20000-1000'
# resName <- 'n20000-p1000-top1000-o0.05-vi0-glmnet 1-5'
# outputDir <- 'D:/projects/robustness/offset-test-10000-200'
# resName <- 'n10000-p200-top1000-o0.1-vi0-glmnet 1-3'

legendLabels <- c('short', 'warm short', 'long', 'warm long')
plotResults(outputDir, resName, legendLabels = legendLabels)
plotRunTime(outputDir, resName, legendLabels = legendLabels)
