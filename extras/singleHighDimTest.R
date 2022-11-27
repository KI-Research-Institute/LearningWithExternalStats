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


n = 3e4  # number of samples in train and tests sets
p = 500  # number of features
binary = T # type of covariates

# Simulation model
outcomeOffset = -log(n/1000)  # offset of the outcome logistic model, determines outcome prevalence
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
    outputDir = getwd(),
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

res <- estimateExternalPerformanceFromStatistics(
  internalData = internalData,
  externalStats = muExt,
  externalEstimatorSettings = externalEstimatorSettings)

cat('Training time:\n')
print(trainingTime)
print(internalAUC, digits=3)
print(extAuc, digits=3)


showResults <- c(
  'n',
  'n outcome',
  'Max weight',
  'chi2 to uniform',
  'kl',
  'Max Weighted SMD',
  'AUROC', '2.5 percentile AUROC', '97.5 percentile AUROC',
  'n repetitions',
  'Brier score',
  'Global calibration mean prediction',
  'Global calibration observed risk')
estimationView <- res$estimation[showResults, , drop = F]
estimationView[, 'value'] <- apply(estimationView, 1, function(x) {sprintf('%.3g', x)})
print(estimationView)
