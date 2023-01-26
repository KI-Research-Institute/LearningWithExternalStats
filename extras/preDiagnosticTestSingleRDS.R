library(glue)
library(LearningWithExternalStats)
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
library(glmnet)
source('./mlWrappers/wglmnet.R')
source('./mlWrappers/wrappedml.R')

n <- 3e5
p <- 500
outputDir <- glue('D:/projects/robustness/offset-test-{n}-{p}')
testName <- 'data_n3e+05-p500-top1000-o0.004-vi0.5_glmnet_10'

dataFileName <- file.path(outputDir, glue('{testName}.rds'))
data <- readRDS(file=dataFileName)
trainingTime <- data$trainingTime
d <- data$d
model1 <- data$model1
xFeatures <- data$xFeatures

vars1 <- wimportant(model1)
cat('Number of selected variables' ,length(vars1), '\n')

# Predict the label probabilities in the internal test set
internalX <- sapply(d$internalTest[xFeatures], as.numeric)
# pInternal <- wpredict(model1, internalX)
# internalAUC <- auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T))

internalK <- d$internalTest[, c(vars1, 'Y')]
externalK <- d$externalTest[, c(vars1, 'Y')]
z <- computeTable1LikeTransformation(internalK, outcomeBalance=TRUE)
dTransformedExt <- computeTable1LikeTransformation(externalK, outcomeBalance=TRUE)
mu <- colMeans(dTransformedExt)

preDiag <- preDiagnostics(z, mu, maxDiff=0.01)
