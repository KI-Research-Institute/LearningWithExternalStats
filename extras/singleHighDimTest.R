library(LearningWithExternalStats)
library(glue)
library(glmnet)
library(pROC)

source('./mlWrappers/wglmnet.R')
source('./mlWrappers/wxgboost.R')
source('./mlWrappers/wrappedml.R')
source('./simulations/anchorModelSimulator.R')


n = 3e5  # number of samples in train and tests sets
p = 500  # number of features
binary = T # type of covariates
ntop = 500  # number of features used in estimation of external performance

# Simulation model
outcomeOffset = -log(250)  # offset of the outcome logistic model, determines outcome prevalence
sigma_B_X_AH = 0  # degree of porximity assumption violation
sigma_B_Y_X_factor = 4
sigma_B_Y_XA_factor = 4
envOffset = 5
# Estimation model
trainer = wglmnet()  # wXGBoost()

# Reweighing parameters
externalEstimatorSettings <-  createExternalEstimatorSettings(
    divergence = 'entropy',  # entropy, chi2
    lambda = 1e-1,
    minW = 0,
    optimizationMethod = 'dual',  # dual, primal
    nMaxReweight = 5000, # maximum samples in a single round of repeated estimations
    stratifiedSampling = T,
    nRepetitions = 10,
    maxProp = 100,
    outputDir = outputDir,
    maxCores = 3
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
print(sort(colMeans(d$internalTest))[floor((0:10)/10*(testParams$p-1)+1)])
trainingStartTime <- Sys.time()
cat('Training using ', testParams$trainer$name, '\n')
model1 <- wfit(testParams$trainer, sapply(d$internalTrain[xFeatures], as.numeric), d$internalTrain[['Y']])
trainingTime <- Sys.time() - trainingStartTime


vars1 <- wimportant(model1)
cat('Number of selected variables' ,length(vars1), '\n')

# Predict the label probabilities in the internal test set
internalX <- sapply(d$internalTest[vars1], as.numeric)
# pInternal <- predict(model1, internalX, type = "response", s = "lambda.1se")[,1]
pInternal <- wpredict(model1, internalX)
internalAUC <- auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T))

# External results
xExternal <- sapply(d$externalTest[vars1], as.numeric)
# pExternal <- predict(model1, xExternal, type = "response", s = "lambda.1se")[,1]
pExternal <- wpredict(model1, xExternal)
extAuc <- auc(roc(d$externalTest[['Y']], pExternal, quiet = TRUE, direction='<'))

observedRisk <- mean(d$externalTest[['Y']])
predictedRisk <- mean(pExternal)
cat('------- Observerved Risk', observedRisk, '----------\n')
cat('------- Predicted   Risk', predictedRisk, '----------\n')

internalY <- d$internalTest[['Y']]
results <- list(
  'n' = length(internalY),
  'n outcome' = sum(internalY),
  'Internal AUC' = internalAUC,
  'External AUC' = extAuc,
  'Training Time' = trainingTime,
  'n Reweight Vars' = length(vars1)
)

internalK <- d$internalTest[, c(vars1, 'Y')]
externalK <- d$externalTest[, c(vars1, 'Y')]
dTransformedInt <- computeTable1LikeTransformation(internalK, outcomeBalance=TRUE)
dTransformedExt <- computeTable1LikeTransformation(externalK, outcomeBalance=TRUE)
muExt <- colMeans(dTransformedExt)
internalData <- list(z=dTransformedInt, p = pInternal, y = internalK[['Y']])
estimationParams <- testParams$estimationParams[[i]]

res <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  externalEstimatorSettings = estimationParams)
