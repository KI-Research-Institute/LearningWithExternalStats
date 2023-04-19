library(glue)
library(LearningWithExternalStats)

testSimpleSimulation <- function(probsInt, probsExt, n, maxProp=100) {
  z <- sapply(probsInt, function(p) {rbinom(n, 1, p)})
  colnames(z) <- paste('X', 1:length(probsInt), sep='')
  # cat('mean z', colMeans(z), '\n')

  zExt <- sapply(probsExt, function(p) {rbinom(n, 1, p)})
  colnames(zExt) <- colnames(z)
  muExt <- colMeans(zExt)
  # cat('mu Ext', muExt, '\n')

  r <- preDiagnostics(z, muExt, 100)

  cat('n represented features:', length(r$representedFeatures), '\n')
  cat('Highly skewed binary:', r$highlySkewedBinary, '\n')
  cat('n =', sum(r$zidx), '\n')
  cat('Status = ', r$status, '\n\n')
  return(r)
}
# TODO Test unary and binary variable check under different scenarios

n <- 3e5
maxProp <- 100

probsInt <- c(0.5, 0.2, 0.1, 0.8,   1, 0.001, rbeta(500, 1, 10))
probsExt <- c(0.5, 0.3, 0,     1, 0.4, 0.2, rbeta(500, 1, 10))
# r1 <- testSimpleSimulation(probsInt, probsExt, n, maxProp)

probsInt <- c(0.5, 0.2, 0.1,  0.8,    1, 0.01)
probsExt <- c(0.5, 0.3, 0,      1,    1, 0.2)
# r2 <- testSimpleSimulation(probsInt, probsExt, n, maxProp)
p <- 100
probsInt <- c(0.5, 0.2, 0.1,  0.8,  0.9, 0.01, rbeta(p, 1, 10))
probsExt <- c(0.5, 0.3, 0,      1,    1, 0.2, rbeta(p, 1, 10))
r3 <- testSimpleSimulation(probsInt, probsExt, n, maxProp)
