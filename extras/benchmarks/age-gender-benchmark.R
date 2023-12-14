#'
#' 1. Read Table 1s in the format of March 23
#' 2. Convert to simple tabular format
#' 3. Complete missing rubrics
#' 4. Check that they are the most common
#' 5. Simulate with/without interatction
#' 6. Emulate the pipeline
#'
#'
#'
#' TODO
#' 1. Make a streamlined pipeline
#' 2. Add number of times there is an estimation
#'
rm(list = ls())

library(pROC)
library(WeightedROC)
library(LearningWithExternalStats)
library(stringr)
library(ggplot2)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
incDir <- '..'
workDir <- 'D:/projects/robustness/age-gender'
# workDir <- 'C:/localdev/projects/robustness/age-gender'

source(file.path(incDir, 'mlWrappers', 'wglmnet.R'))
source(file.path(incDir, 'mlWrappers', 'wxgboost.R'))
source(file.path(incDir, 'mlWrappers', 'wrappedml.R'))
source(file.path(incDir, 'optimize', 'kerasWeightOptimizer.R'))


source('./age-gender-utils.R')

allDbNames <- c('ccae', 'mdcd', 'mdcr', 'optum.claims', 'optum.ehr')

nMaxTrainSample <- 5000  # 50000
nSimulationRepetitions <- 5  # number of repetitions per configuration (5)
nEstimationRepetitions <- 10  # number of repetitions in estimation
maxCores <- 10

outputDir <- file.path(workDir, 'output')
reweightEg500 <- seTunedWeightOptimizer(outcomeCol = 'Y', outputDir=outputDir, nIter = 500)
reweightEg2000 <- seTunedWeightOptimizer(outcomeCol = 'Y', outputDir=outputDir, nIter = 2000)
reweightAdamId <- kerasWeightOptimizer(optimizer = 'adam', nIter = 10000, parameterization = 'id')
reweightAdamLog <- kerasWeightOptimizer(optimizer = 'adam', nIter = 10000, parameterization = 'log')
reweightSGDId <- kerasWeightOptimizer(optimizer = 'sgd', nIter = 10000, parameterization = 'id')

reweightConfigs <- list(
  list('Table 1', reweightEg500, 'EG C 500'),
  list('Table 1', reweightEg2000, 'EG C 2000'),
  list('Flat',  reweightEg2000, 'EG N'),
  # list('New flat',  reweightAdamId, 'Adam Id N'),
  # list('New flat',  reweightSGDId, 'SGD Id N'),
  # list('New flat',  reweightKerasLog, 'KerasLog N'),
  list('Interaction v0',  reweightEg2000, 'EG Iv0'),
  list('Interaction',  reweightEg2000, 'EG I')
  # list('Interaction',  reweightAdamId, 'Adam Id I'),
  # list('Interaction',  reweightSGDId, 'SGD Id I')
  # list('Interaction',  reweightKerasLog, 'KerasLog I')
)

runTests <- T

skewNames <- c('', 'skew 4', 'skew 3', 'skew 1')
# skewNames <- c('', 'cskew 4')

analysisIdxs <- c(1, 3, 5, 7, 9)  #
# analysisIdxs <- c(1, 3)


for (trainer in list(wglmnet())) {
  for (prefix in c('', 'clinical64 ')) {  # , 'clinical64-'
    # Main loop
    for (analysisIdx in analysisIdxs) {
      if (runTests) {
        for (skewName in skewNames) {
        if (nchar(skewName)>0) {
          skew <- read.csv(file.path(workDir, 'input', glue('{skewName}.csv')), row.names = 'X')
        } else {
          skew <- NULL
        }
        table1s <- read.csv(
          file.path(workDir, 'output', 'short table 1', glue('{prefix} Table 1 - {analysisIdx}.csv')), row.names = 'X')
        tag <- getSkewSpecificTag(prefix, trainer, analysisIdx, nMaxTrainSample, nSimulationRepetitions, skewName)
        cat(tag, '\n')
        rdsFileName <- file.path(workDir, 'output', glue('{tag}-results.RDS'))
        allResults <- vector(mode = 'list', length = nSimulationRepetitions)
        nDbs <- length(allDbNames)
        for (b in 1:nSimulationRepetitions)
          for (intDbIdx in 1:nDbs) {
            intDbName <- allDbNames[intDbIdx]
            for (extDbIdx in 1:nDbs) { # nDbs
              extDbName <- allDbNames[extDbIdx]
              if (intDbIdx != extDbIdx) {
                exprimentName <- glue(
                  '{prefix}{allDbNames[intDbIdx]}-Analysis_{analysisIdx}-{allDbNames[extDbIdx]}-{b}')
                cat('\n', exprimentName, '\n')
                allResults[[exprimentName]] <- testSingleAgeGenderOutcome(
                  table1s, intDbName, extDbName, reweightConfigs, nMaxTrainSample = nMaxTrainSample, skew = skew,
                  trainer = trainer, nRepetitions=nEstimationRepetitions, maxCores=maxCores, outputDir=outputDir)
              }
            }
          }
        saveRDS(allResults, file.path(rdsFileName))
        rSummary <- arrangeAgeGenderEmulationResults(allResults, reweightConfigs)
        p <- rSummary %>%
          ggplot(aes(x=type, y=diff, color = type)) + # group=type,
          geom_boxplot()
        ggsave(file.path(workDir, 'output', glue('{tag} ext-eval.png')))
        }
      }
      plotTable1Results(workDir, prefix, skewNames, nMaxTrainSample, nSimulationRepetitions, analysisIdx)
    }
  }
}
