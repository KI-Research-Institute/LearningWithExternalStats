#' Get default Model Hyper-parameters
#'
#' @param p covariates dimension
#' @param sigma_B_X_AH the degree of proximity assumption violation (0 none, 1 high)
#'
#' @return an object of AnchorHyperParams
#'
getDefaultModelHyperParams <- function(
    p=500, sigma_B_X_AH=0, outcomeOffset=0, sigma_B_Y_X_factor=4, sigma_B_Y_XA_factor=4) {
  hyperParams = list(
    p = p,
    n_binary = 0,
    sigma_M_H = 0.2,
    sigma_M_X = 0.2,  # std of M_X, coefficient of the influence of A on X
    sigma_M_Y = 0.2,

    sigma_B_X_H = 1,
    sigma_B_Y_H = 1,

    sigma_B_X_AH = sigma_B_X_AH, # induce different correlation per environment
    sigma_B_Y_X = sigma_B_Y_X_factor/sqrt(p),
    sigma_B_Y_XA = sigma_B_Y_XA_factor/sqrt(p),

    outcomeOffset = outcomeOffset,
    shape1 = 1,
    shape2 = 10
  )
  class(hyperParams) <- 'AnchorHyperParams'
  return(hyperParams)
}


#' Get model parameters
#'
#' @param hyperParams hyper parameters of class 'AnchorHyperParams'
#'
#' @return an object of AnchorModelParams
#'
modelParams <- function(hyperParams) {
  # TODO merge with anchor model
  params <- list(

    q = 1,  # dim(H)
    p = hyperParams$p, # dim(X)

    A = NaN,

    CeX = diag(hyperParams$p),  # Covariance matrix of epsilon_X
    M_H = rnorm(1, 0, hyperParams$sigma_M_H),

    M_Y = rnorm(1, 0, hyperParams$sigma_M_Y),

    M_X = rnorm(hyperParams$p, 0, hyperParams$sigma_M_X),  # coefficients of influence of A on X
    B_X_H = rnorm(hyperParams$p, 0, hyperParams$sigma_B_X_H),
    B_X_AH = rnorm(hyperParams$p, 0, hyperParams$sigma_B_X_AH),  # induce different correlation per environment

    B_Y_X = rnorm(hyperParams$p, 0, hyperParams$sigma_B_Y_X),
    B_Y_H = rnorm(1, 0, hyperParams$sigma_B_Y_H),
    B_Y_XA = rnorm(hyperParams$p, 0, hyperParams$sigma_B_Y_XA),

    outcomeOffset = hyperParams$outcomeOffset,

    n_binary = hyperParams$n_binary,
    # in case of binary features frequencies are taken from a beta distribution
    freqs = rbeta(hyperParams$p, hyperParams$shape1, hyperParams$shape2)
  )
  class(params) <- 'AnchorModelParams'
  return(params)
}

#' Get model parameters
#'
#' @param m model parameters, an object of AnchorModelParams
#' @param n number of samples
#'
#' @return a data frame
#'
sampleModel <- function(m, n)
{
  xFeatures <- glue(('X{1:m$p}'))
  p <- m$p
  q <- m$q

  d <- data.frame(matrix(ncol = p+1, nrow=n))
  colnames(d)[p+1] <- 'Y'
  colnames(d)[1:p] <- xFeatures

  H <- m$A *  m$M_H + matrix( rnorm(n*q,mean=0,sd=1), n, q)

  e_X = matrix( rnorm(n*p,mean=0,sd=1), n, p) %*% m$CeX
  X <- m$A * rep(1, n) %*% t(m$M_X) + H %*% m$B_X_H + m$A * H %*% m$B_X_AH  + e_X  #
  if (m$n_binary>0) {
    n_binary <- min(m$n_binary, p)
    for (i in (p-n_binary+1):p) {
      probHiddenXi <- pnorm(X[, i], mean(X[, i]) ,sd(X[, i]))
      X[, i] <- as.numeric(probHiddenXi < m$freqs[i])
    }
  }
  colnames(X) <- xFeatures

  d[ ,xFeatures] <- X
  logit <- (X %*% m$B_Y_X + m$A * X %*%m$B_Y_XA + m$B_Y_H* H + m$M_Y * m$A + m$outcomeOffset)
  pY <- 1/(1+exp(-logit))
  d[ ,'Y'] <- rbinom(n, size = 1, prob = pY)

  return (d)
}


generateSimulatedData <- function(testParams) {
  hyperParams = getDefaultModelHyperParams(
    p=testParams$p, outcomeOffset = testParams$outcomeOffset, sigma_B_X_AH = testParams$sigma_B_X_AH,
    sigma_B_Y_X_factor = testParams$sigma_B_Y_X_factor, sigma_B_Y_XA_factor = testParams$sigma_B_Y_XA_factor)
  hyperParams$n_binary <- testParams$n_binary
  params <- modelParams(hyperParams)
  params$A <- 0
  internalTrain <- sampleModel(params, testParams$n)
  internalTest <- sampleModel(params, testParams$n)
  params$A <- testParams$envOffset
  externalTest <- sampleModel(params, testParams$n)
  d <- list(internalTrain = internalTrain, internalTest = internalTest, externalTest = externalTest)
  cat(glue("n outcome {sum(internalTest['Y'])}/{nrow(internalTest)}"),"\n")
  return(d)
}
