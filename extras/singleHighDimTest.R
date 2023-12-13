rm(list = ls())
gc()
if (!is.null(dev.list()["RStudioGD"]))
  dev.off(dev.list()["RStudioGD"])

library(LearningWithExternalStats)
library(glue)
library(glmnet)
library(pROC)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)

source('./mlWrappers/wglmnet.R')
source('./mlWrappers/wxgboost.R')
source('./mlWrappers/wrappedml.R')
source('./simulations/anchorModelSimulator.R')
source('../../explore/FastBootstrap.R')


n = 1000  # number of samples in train and tests sets
p = 10 # number of features
binary = T # type of covariates
lambda = 1e-1

# Optimization algorithms
# wOptimizers <- vector(mode = 'list', length = 2)

minEpsilon0 <- 1e-4
maxEpsilon0 <- 1e-1
nEpsilon0 <- 10


wOptimizer <- cyclicOptimizer(
  minSd = 1e-4,
  nIter = 40,
  p = 1e-1,
  outputDir = 'D:/projects/robustness/high-dim-small'
)



wOptimizer <- cvxWeightOptimizer()



minAlpha <- 0.1
maxAlpha <- 10
nAlphas <- 5
wOptimizer <- seTunedWeightOptimizer(
  outcomeCol = 'Y',
  nIter = 1000,
  alphas = minAlpha*exp(log(maxAlpha/minAlpha)*(0:nAlphas)/nAlphas),
  improveTh = 1,
  maxErr = 1e-3,
  approxUpdate = F,
  outputDir = 'D:/projects/robustness/high-dim-small'
)


wOptimizer <- sgdTunedWeightOptimizer(
  epsilon0s = minEpsilon0*exp(log(maxEpsilon0/minEpsilon0)*(0:nEpsilon0)/nEpsilon0),
  polyBeta = 1,
  batchSize = 100,
  nIter = 1000,
  nProbe = 50,
  normalizeNu = T,
  outputDir = 'D:/projects/robustness/high-dim-small',
  averaging = T
)




# Simulation model
outcomeOffset = -log(4)  # offset of the outcome logistic model, determines outcome prevalence
sigma_B_X_AH = 0  # degree of porximity assumption violation
sigma_B_Y_X_factor = 4
sigma_B_Y_XA_factor = 4
envOffset = 5
# Estimation model
trainer = wglmnet()  # wXGBoost()

# Reweighing parameters
externalEstimatorSettings <-  createExternalEstimatorSettings(
  reweightAlgorithm = wOptimizer,
  nMaxReweight = n, # maximum samples in a single round of repeated estimations
  stratifiedSampling = T,
  nRepetitions = 200,
  maxProp = 100,
  outputDir = getwd(),
  maxCores = 15
)

# Simulate data
hyperParams = getDefaultModelHyperParams(
  p=p, outcomeOffset = outcomeOffset, sigma_B_X_AH = sigma_B_X_AH,
  sigma_B_Y_X_factor = sigma_B_Y_X_factor, sigma_B_Y_XA_factor = sigma_B_Y_XA_factor)
hyperParams$binary <- binary
params <- modelParams(hyperParams)
params$A <- 0
internalTrain <- sampleModel(params, n)
internalTest <- sampleModel(params, n)
params$A <- envOffset
externalTest <- sampleModel(params, n)
d <- list(internalTrain = internalTrain, internalTest = internalTest, externalTest = externalTest)
cat(glue("n outcome {sum(internalTest['Y'])}/{nrow(internalTest)}"),"\n")


xFeatures <- colnames(d$internalTest)[1:(ncol(d$internalTest)-1)]
cat('Generated simulated data, with the following means:\n')
print(sort(colMeans(d$internalTest))[floor((0:10)/10*(p-1)+1)])
trainingStartTime <- Sys.time()
cat('Training using ', trainer$name, '\n')
model1 <- wfit(trainer, sapply(d$internalTrain[xFeatures], as.numeric), d$internalTrain[['Y']])
trainingTime <- Sys.time() - trainingStartTime


vars1 <- wimportant(model1)
cat('Number of selected variables' ,length(vars1), '\n')

# Predict the label probabilities in the internal test set
internalX <- sapply(d$internalTest[xFeatures], as.numeric)
# pInternal <- predict(model1, internalX, type = "response", s = "lambda.1se")[,1]
pInternal <- wpredict(model1, internalX)
internalAUC <- auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T))

# External results
xExternal <- sapply(d$externalTest[xFeatures], as.numeric)
# pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
pExternal <- wpredict(model1, xExternal)
extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))

observedRisk <- mean(d$externalTest[['Y']])
predictedRisk <- mean(pExternal)

internalY <- d$internalTest[['Y']]

internalK <- d$internalTest[, c(vars1, 'Y')]
externalK <- d$externalTest[, c(vars1, 'Y')]
dTransformedInt <- computeTable1LikeTransformation(internalK, outcomeBalance=TRUE)
dTransformedExt <- computeTable1LikeTransformation(externalK, outcomeBalance=TRUE)
muExt <- colMeans(dTransformedExt)
internalData <- list(z=dTransformedInt, p = pInternal, y = internalK[['Y']])


result <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  externalEstimatorSettings = externalEstimatorSettings)

comparison <- list(
  extAuc=extAuc,
  internalAUC=internalAUC,
  estAUC=result$estimation['AUROC', ])

cat('\nSummary:\n')
print(unlist(comparison), digits=3)
cat('Estimation time ')
print(result$estimationTime)

print(result$estimation)



