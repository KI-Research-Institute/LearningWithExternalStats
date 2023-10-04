library(pROC)


test_that("real AUCs are consistent", {
  # Load package data
  d <- LearningWithExternalStats::anchorData1
  model1 <- LearningWithExternalStats::anchorLR1
  pInternal <- predictGLM(d$internalTest, model1)
  internalAUC <- as.numeric(auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T)))
  expect_equal(internalAUC, 0.93, tolerance = 0.02)

  pExternal <- predictGLM(d$externalTest, model1)
  externalAUC <- as.numeric(auc(roc(d$externalTest[['Y']], pExternal, direction = "<", quiet = T)))
  expect_equal(externalAUC, 0.72, tolerance = 0.05)

})


test_that("estimated AUC is in range", {

  replaceSimpleLoggerThreshold('WARN')

  d <- LearningWithExternalStats::anchorData1
  model1 <- LearningWithExternalStats::anchorLR1
  pInternal <- predictGLM(d$internalTest, model1)
  dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
  dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
  muExt <- colMeans(dTransformedExt)

  pExternal <- predictGLM(d$externalTest, model1)
  externalAUC <- as.numeric(auc(roc(d$externalTest[['Y']], pExternal, direction = "<", quiet = T)))

  wOptimizer = seTunedWeightOptimizer(maxSuccessMSE=5e-4)
  estimatedAUC <- estimateExternalAUCWithInternalData(dTransformedInt, pInternal, muExt, wOptimizer)
  expect_equal(estimatedAUC, externalAUC, tolerance = 0.07)

})


# Robustness to exclusions


test_that("algorithm detects missing stats", {

  replaceSimpleLoggerThreshold('WARN')

  excludeVars <- c('X3:y1', 'X3:y0')
  d <- LearningWithExternalStats::anchorData1
  model1 <- LearningWithExternalStats::anchorLR1
  pInternal <- predictGLM(d$internalTest, model1)
  dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
  dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
  muExt <- colMeans(dTransformedExt)
  reducedMu <- muExt[!names(muExt) %in% excludeVars]
  expect_equal(length(muExt)-2, length(reducedMu))

  pExternal <- predictGLM(d$externalTest, model1)
  externalAUC <- as.numeric(auc(roc(d$externalTest[['Y']], pExternal, direction = "<", quiet = T)))

  wOptimizer <- seTunedWeightOptimizer()
  estimatedAUC <- estimateExternalAUCWithInternalData(dTransformedInt, pInternal, reducedMu, wOptimizer)
  expect_null(estimatedAUC)

})


test_that("algorithm detects missing internal features", {

  replaceSimpleLoggerThreshold('WARN')

  excludeVars <- c('X3:y1', 'X3:y0')
  d <- LearningWithExternalStats::anchorData1
  model1 <- LearningWithExternalStats::anchorLR1
  pInternal <- predictGLM(d$internalTest, model1)
  dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
  dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
  muExt <- colMeans(dTransformedExt)
  zReduced <- dTransformedInt[, !colnames(dTransformedInt) %in% excludeVars]

  pExternal <- predictGLM(d$externalTest, model1)
  externalAUC <- as.numeric(auc(roc(d$externalTest[['Y']], pExternal, direction = "<", quiet = T)))

  wOptimizer <- seTunedWeightOptimizer()
  estimatedAUC <- estimateExternalAUCWithInternalData(zReduced, pInternal, muExt, wOptimizer)
  expect_null(estimatedAUC)
})


## Binary data tests ###############################################################################################


test_that("real AUCs are consistent with binary data", {

  replaceSimpleLoggerThreshold('WARN')

  # Load package data
  d <- LearningWithExternalStats::binaryAnchorData1
  model1 <- LearningWithExternalStats::binaryAnchorLR1
  pInternal <- predictGLM(d$internalTest, model1)
  internalAUC <- as.numeric(auc(roc(d$internalTest[['Y']], pInternal, direction = "<", quiet = T)))
  cat(internalAUC, '\n')
  expect_equal(internalAUC, 0.72, tolerance = 0.02)

  pExternal <- predictGLM(d$externalTest, model1)
  externalAUC <- as.numeric(auc(roc(d$externalTest[['Y']], pExternal, direction = "<", quiet = T)))
  expect_equal(externalAUC, 0.65, tolerance = 0.02)

})


test_that("estimated AUC is in range with binary data", {

  replaceSimpleLoggerThreshold('WARN')

  d <- LearningWithExternalStats::binaryAnchorData1
  model1 <- LearningWithExternalStats::binaryAnchorLR1
  pInternal <- predictGLM(d$internalTest, model1)
  dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
  dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
  muExt <- colMeans(dTransformedExt)

  pExternal <- predictGLM(d$externalTest, model1)
  externalAUC <- as.numeric(auc(roc(d$externalTest[['Y']], pExternal, direction = "<", quiet = T)))

  wOptimizer = seTunedWeightOptimizer()
  estimatedAUC <- estimateExternalAUCWithInternalData(dTransformedInt, pInternal, muExt, wOptimizer)
  expect_equal(estimatedAUC, externalAUC, tolerance = 0.05)

})


test_that("package caputures the status", {

  replaceSimpleLoggerThreshold('WARN')

  d <- LearningWithExternalStats::binaryAnchorData1
  model1 <- LearningWithExternalStats::binaryAnchorLR1
  pInternal <- predictGLM(d$internalTest, model1)
  dTransformedInt <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
  dTransformedExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
  muExt <- colMeans(dTransformedExt)

  pExternal <- predictGLM(d$externalTest, model1)
  externalAUC <- as.numeric(auc(roc(d$externalTest[['Y']], pExternal, direction = "<", quiet = T)))

  wOptimizer = seTunedWeightOptimizer(maxSuccessMSE=5e-4)
  estimation1 <- estimatePerformanceWithInternalData(dTransformedInt, pInternal, muExt, wOptimizer)
  expect_equal(estimation1$status=='Success', T)

  wOptimizer = seTunedWeightOptimizer(nTuneIter = 2, nIter = 5, maxSuccessMSE=0.001)
  estimation2 <- estimatePerformanceWithInternalData(dTransformedInt, pInternal, muExt, wOptimizer)
  expect_equal(estimation2$status=='Success', F)

})




