test_that("binary data passes", {
  d <- LearningWithExternalStats::binaryAnchorData1
  z <- d$internalTest
  mu <- colMeans(d$externalTest)
  preDiag <- preDiagnostics(z, mu, maxDiff=0.01)
  expect_equal(preDiag$status == 'Success', T)
})


test_that("diagnostics detects NA in Z and mu", {
  d <- LearningWithExternalStats::binaryAnchorData1
  z <- d$internalTest
  mu <- colMeans(d$externalTest)
  muNA <- mu
  muNA[3] <- NA
  preDiag <- preDiagnostics(z, muNA, maxDiff=0.01)
  expect_equal(preDiag$status == 'Success', F)
  z[10,11] <- NA
  preDiag <- preDiagnostics(z, mu, maxDiff=0.01)
  expect_equal(preDiag$status == 'Success', F)
})


test_that("skewed binary features are detected", {
  d <- LearningWithExternalStats::binaryAnchorData1
  z <- d$internalTest
  mu <- colMeans(d$externalTest)
  n <- nrow(z)
  z[, 1] <- rbinom(n, 1, 5/n)
  mu[1] <- 0.5
  preDiag <- preDiagnostics(z, mu, maxDiff=0.01)
  expect_equal(preDiag$status == 'Success', F)
})


test_that("skewed unary features are detected", {
  d <- LearningWithExternalStats::binaryAnchorData1
  z <- computeTable1LikeTransformation(d$internalTest, outcomeBalance=TRUE)
  zExt <- computeTable1LikeTransformation(d$externalTest, outcomeBalance=TRUE)
  mu <- colMeans(zExt)
  preDiag <- preDiagnostics(z, mu, maxDiff=0.01)
  expect_equal(preDiag$status == 'Success', T)
  mu['X1_Table1T_times_y1'] <- 0.011
  preDiag <- preDiagnostics(z, mu, maxDiff=0.01)
  expect_equal(preDiag$status == 'Success', F)

})
