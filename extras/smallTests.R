rm(list = ls())

library(LearningWithExternalStats)
library(glue)
library(pROC)
library(glmnet)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./optimize/kerasWeightOptimizer.R')


d <- LearningWithExternalStats::binaryAnchorData1
model1 <- LearningWithExternalStats::binaryAnchorLR1


xFeatures <- colnames(d$externalTest)[1:(ncol(d$externalTest)-1)]
internalX <- sapply(d$internalTest[xFeatures], as.numeric)
pInternal <- predict(model1, internalX, type = "response", s = "lambda.1se")[,1]
internalAUC <- auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T))
cat(glue('\nInternal AUC = {format(internalAUC, digits=3)}'), '\n')

xExternal <- sapply(d$externalTest[xFeatures], as.numeric)
pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))
cat(glue('\nExternal AUC = {format(extAuc, digits=3)}'), '\n')


# Table1Transform <- reweightTransfrom$new(outcomeCol = 'Y', interactionVars =  glue('X{1:2}'))  # , interactionVars = glue('X{1:2}')
# f <- Table1Transform$getFormula(d$internalTest)

# dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
# dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)

dTransformedInt <-transformClassifierData(
  d$internalTest,
  transformType = 'Interaction',
  interactionVars = glue('X{1:2}'),
  outcomeCol = 'Y'
)
dTransformedExt <- transformClassifierData(
  d$externalTest,
  transformType = 'Interaction',
  interactionVars = glue('X{1:2}'),
  outcomeCol = 'Y'
)

muExt <- colMeans(dTransformedExt)


externalEstimatorSettings <- createExternalEstimatorSettings(
  reweightAlgorithm =  kerasWeightOptimizer(optimizer = 'adam', nIter = 10000, parameterization = 'id'),
  # reweightAlgorithm = seTunedWeightOptimizer(nIter = 5000, maxSuccessMSE = 1e-4), # cvxWeightOptimizer(),  nIter = 10000
  nMaxReweight = 10000,
  nRepetitions = 10,
  maxCores = 1,
  maxWSMD = 0.1
)

internalData <- list(z=dTransformedInt, p = pInternal, y = d$internalTest[['Y']])

estimatedLRResults1 <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  externalEstimatorSettings = externalEstimatorSettings
)

estimationLRView <- estimatedLRResults1$estimation
if (!is.null(estimationLRView)) {
  estimationLRView[, 'value'] <- apply(estimationLRView, 1, function(x) {sprintf('%.3g', x)})
  # print(estimationLRView)
} else
{
  cat('Empty estimation results\n')
}


fields <- c(getPreDiagnosticsFieldNames(), getWeightingResultsFieldNames(), getEstimationFieldNames(),
            'Internal AUC', 'External AUC')
summary <- data.frame(row.names = fields)

expName <- 'Exp.1'

summary['External AUC', expName] <- extAuc
summary['Internal AUC', expName] <- internalAUC
summary[rownames(estimatedLRResults1$results), expName] <- estimatedLRResults1$results
# write.csv(estimatedLRResults1$results, 'results1.csv')
cat('Estimated results', estimatedLRResults1$estimation['AUROC', 'value'], '\n')


# Test dropping a feature
expName <- 'Exp.2'
cat('mu length =', length(muExt), '\n')
excludeVars <- c('I(Y):X3', 'X3:I(1 - Y)')
reducedMu <- muExt[!names(muExt) %in% excludeVars]

cat('Reduced mu length =', length(reducedMu), '\n')
estimatedLRResults2 <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = reducedMu,
  externalEstimatorSettings = externalEstimatorSettings
)
cat('Estimated results with and excluded var from stats', estimatedLRResults2$estimation['AUROC', 'value'], '\n')
summary[rownames(estimatedLRResults2$results), expName] <- estimatedLRResults2$results
write.csv(estimatedLRResults2$results, 'results2.csv')

expName <- 'Exp.3'
internalData$z <- internalData$z[, !colnames(internalData$z) %in% excludeVars]
cat('Reduced z width =', ncol(internalData$z), '\n')
estimatedLRResults3 <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  externalEstimatorSettings = externalEstimatorSettings
)
cat('Estimated results with and excluded var from data', estimatedLRResults3$estimation['AUROC', 'value'], '\n')
summary[rownames(estimatedLRResults3$results), expName] <- estimatedLRResults3$results
write.csv(estimatedLRResults3$results, 'results3.csv')
