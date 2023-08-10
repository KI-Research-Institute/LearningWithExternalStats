rm(list = ls())

library(LearningWithExternalStats)
library(glue)
library(pROC)
library(glmnet)

d <- LearningWithExternalStats::binaryAnchorData1
model1 <- LearningWithExternalStats::binaryAnchorLR1


xFeatures <- colnames(d$externalTest)[1:(ncol(d$externalTest)-1)]
internalX <- sapply(d$internalTest[xFeatures], as.numeric)
pInternal <- predict(model1, internalX, type = "response", s = "lambda.1se")[,1]
internalAUC <- auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T))
cat(glue('\nInternal AUC = {format(internalAUC, digits=3)}'), '\n')

dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
muExt <- colMeans(dTransformedExt)


externalEstimatorSettings <- createExternalEstimatorSettings(
  reweightAlgorithm = seTunedWeightOptimizer(maxSuccessMSE = 1e-4), # cvxWeightOptimizer(),
  nMaxReweight = 10000,
  nRepetitions = 10,
  maxCores = 3,
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
  print(estimationLRView)
} else
{
  cat('Empty estimation results\n')
}


fields <- c(getPreDiagnosticsFieldNames(), getWeightingResultsFieldNames(), getEstimationFieldNames(),
            'Internal AUC', 'External AUC')
summary <- data.frame(row.names = fields)

expName <- 'Exp.1'
xExternal <- sapply(d$externalTest[xFeatures], as.numeric)
pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))
cat(glue('\nExternal AUC = {format(extAuc, digits=3)}'), '\n')

summary['External AUC', expName] <- extAuc
summary['Internal AUC', expName] <- internalAUC
summary[rownames(estimatedLRResults1$results), expName] <- estimatedLRResults1$results
write.csv(estimatedLRResults1$results, 'results1.csv')


# Test dropping a feature
expName <- 'Exp.2'
cat('mu length =', length(muExt), '\n')
excludeVars <- c('X3_Table1T_times_y1', 'X3_Table1T_times_y0')
reducedMu <- muExt[!names(muExt) %in% excludeVars]

cat('Reduced mu length =', length(reducedMu), '\n')
estimatedLRResults2 <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = reducedMu,
  externalEstimatorSettings = externalEstimatorSettings
)
cat(estimatedLRResults2$estimation['AUROC', 'value'], '\n')
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
cat(estimatedLRResults3$estimation['AUROC', 'value'], '\n')
summary[rownames(estimatedLRResults3$results), expName] <- estimatedLRResults3$results
write.csv(estimatedLRResults3$results, 'results3.csv')
