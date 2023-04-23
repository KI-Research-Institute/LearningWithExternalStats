rm(list=ls())

library(glue)
library(LearningWithExternalStats)

testSimpleSimulation <- function(probsInt, probsExt, n) {
  z <- sapply(probsInt, function(p) {rbinom(n, 1, p)})
  colnames(z) <- paste('X', 1:length(probsInt), sep='')

  zExt <- sapply(probsExt, function(p) {rbinom(n, 1, p)})
  colnames(zExt) <- colnames(z)
  muExt <- colMeans(zExt)

  r <- preDiagnostics(z, muExt, 0.01)

  cat('n represented features:', length(r$representedFeatures), '\n')
  cat('Highly skewed binary:', r$highlySkewedBinary, '\n')
  cat('n =', sum(r$zidx), '\n')
  cat('Status = ', r$status, '\n\n')
  print(r$structuredLog)
  return(r)
}
# TODO Test unary and binary variable check under different scenarios

n <- 1e4

if (T) {
  probsInt <- c(0.5, 0.2, 0.1, 0.8, 0.9, 0.001, rbeta(500, 1, 10))
  probsExt <- c(0.5, 0.3, 0,     1, 0.4, 0.2, rbeta(500, 1, 10))
  r1 <- testSimpleSimulation(probsInt, probsExt, n)
}
if (T) {
  probsInt <- c(0.5, 0.2, 0.1,  0.8,    1, 0.01)
  probsExt <- c(0.5, 0.3, 0,      1,    1, 0.2)
  r2 <- testSimpleSimulation(probsInt, probsExt, n)
}
if (T) {
  probsInt <- c(0.5, 0.2, 0.1, 0.8,   1, 1e-6)
  probsExt <- c(0.5, 0.3, 0,     1, 0.4, 0.2)
  r3 <- testSimpleSimulation(probsInt, probsExt, n)
}
if (T) {
  probsInt <- c(0.5, 0.2, 0.1,  0.8,    1, 0.1)
  probsExt <- c(0.5, 0.3, 0,      1,    1, 0.2)
  r2 <- testSimpleSimulation(probsInt, probsExt, n)
}

#p <- 100
#probsInt <- c(0.5, 0.2, 0.1,  0.8,  0.9, 0.01, rbeta(p, 1, 10))
#probsExt <- c(0.5, 0.3, 0,      1,    1, 0.2, rbeta(p, 1, 10))
#r3 <- testSimpleSimulation(probsInt, probsExt, n)
