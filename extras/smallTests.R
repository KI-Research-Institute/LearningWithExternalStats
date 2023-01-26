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
  reweightAlgorithm = seTunedWeightOptimizer(), # cvxWeightOptimizer(),
  nMaxReweight = 10000,
  nRepetitions = 1,
  maxCores = 1
)

internalData <- list(z=dTransformedInt, p = pInternal, y = d$internalTest[['Y']])

estimatedLRResults <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  externalEstimatorSettings = externalEstimatorSettings
)

showResults <- c(
  'n',
  'n outcome',
  'AUROC',
  'Brier score',
  'Global calibration mean prediction',
  'Global calibration observed risk')
estimationLRView <- estimatedLRResults$estimation[showResults, , drop = F]
estimationLRView[, 'value'] <- apply(estimationLRView, 1, function(x) {sprintf('%.3g', x)})
print(estimationLRView)


xExternal <- sapply(d$externalTest[xFeatures], as.numeric)
pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))
cat(glue('\nExternal AUC = {format(extAuc, digits=3)}'), '\n')


# Test dropping a feature
cat('mu length =', length(muExt), '\n')
excludeVars <- c('X3_Table1T_times_y1', 'X3_Table1T_times_y0')
reducedMu <- muExt[!names(muExt) %in% excludeVars]

cat('Reduced mu length =', length(reducedMu), '\n')
estimatedLRResults <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = reducedMu,
  externalEstimatorSettings = externalEstimatorSettings
)
cat(estimatedLRResults$estimation['AUROC', 'value'], '\n')


internalData$z <- internalData$z[, !colnames(internalData$z) %in% excludeVars]
cat('Reduced z width =', ncol(internalData$z), '\n')
estimatedLRResults <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  externalEstimatorSettings = externalEstimatorSettings
)
cat(estimatedLRResults$estimation['AUROC', 'value'], '\n')
