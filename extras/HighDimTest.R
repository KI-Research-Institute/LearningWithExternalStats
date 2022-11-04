rm(list=ls())
library(LearningWithExternalStats)
library(glue)
library(glmnet)
library(pROC)
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('../data-raw/anchorModel.R')

n <- 50000

# Generate simulated data
hyperParams = getDefaultAnchorHyperParams()
hyperParams$sigma_B_X_AH = 0.5  # Medium shift
params <- anchorParams(hyperParams)
params$A <- 0
internalTrain <- sampleAnchor(params, n)
internalTest <- sampleAnchor(params, n)
params$A <- 1
externalTest <- sampleAnchor(params, n)
d <- list(internalTrain = internalTrain, internalTest = internalTest, externalTest = externalTest)

# Train a model
xFeatures <- colnames(d$internalTest)[1:(ncol(d$internalTest)-1)]
model1 <- cv.glmnet(sapply(d$internalTrain[xFeatures], as.numeric), d$internalTrain[['Y']],
                    family = "binomial", type.measure = "auc", alpha = 0.5)

# Predict the label probabilities in the internal test set
internalX <- sapply(d$internalTest[xFeatures], as.numeric)
pInternal <- predict(model1, internalX, type = "response", s = "lambda.1se")[,1]
internalAUC <- auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T))
cat(glue('\nInternal AUC = {format(internalAUC, digits=3)}'), '\n')

## Estimation using reweighting
dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
muExt <- colMeans(dTransformedExt)
reweightSettings <- createReweightSettings(
  divergence = 'entropy',
  lambda = 1e-2,
  minW = 1e-6,
  optimizationMethod = 'dual'
)
internalData <- list(z=dTransformedInt, p = pInternal, y = d$internalTest[['Y']])
estimatedResults <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  reweightSettings = reweightSettings,
  nboot = 0
)
print(format(estimatedResults$summary, digits=3))

# External results
xExternal <- sapply(d$externalTest[xFeatures], as.numeric)
pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))
cat(glue('\nExternal AUC = {format(extAuc, digits=3)}'), '\n')
