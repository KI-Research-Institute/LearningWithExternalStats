
getDefaultAnchorHyperParams <- function() {
  hyperParams = list(
    p = 10,
    sigma_M_H = 0.2,
    sigma_M_X = 0.2,  # std of M_X, coefficient of the influence of A on X
    sigma_M_Y = 0.2,

    sigma_B_X_H = 1,
    sigma_B_Y_H = 1,

    sigma_B_X_AH = NaN, # induce different correlation per envrionment
    B_Y_X = c(1, 1),
    B_Y_XA = c(-0.8, -0.2)  # May influence overlap in specific classes, difference in classider oreintation

  )
  return(hyperParams)
}


anchorParams <- function(hyperParams) {
  params <- list(

    q = 1,  # dim(H)
    p = hyperParams$p, # dim(X)

    A = NaN,

    CeX = diag(hyperParams$p),  # Covariance matrix of epsilon_X
    M_H = rnorm(1, 0, hyperParams$sigma_M_H),

    M_Y = rnorm(1, 0, hyperParams$sigma_M_Y),

    M_X = rnorm(hyperParams$p, 0, hyperParams$sigma_M_X),  # coefficients of influence of A on X
    B_X_H = rnorm(hyperParams$p, 0, hyperParams$sigma_B_X_H),
    B_X_AH = rnorm(hyperParams$p, 0, hyperParams$sigma_B_X_AH),  # induce different correlation per envrionment

    B_Y_X = c(hyperParams$B_Y_X, rep(0, hyperParams$p-2)),
    B_Y_H = rnorm(1, 0, hyperParams$sigma_B_Y_H),
    B_Y_XA = c(hyperParams$B_Y_XA, rep(0, hyperParams$p-2)) # smaller effect
  )
  return(params)
}

# TODO ALLOW DIFFERENT COVARIANCE PER OUTCOME
# Covaiance matrices

# Generate internal and external data
sampleAnchor <- function(m, n)
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
  colnames(X) <- xFeatures

  d[ ,xFeatures] <- X
  logit <- (X %*% m$B_Y_X + m$A * X %*%m$B_Y_XA + m$B_Y_H* H + m$M_Y * m$A)
  pY <- 1/(1+exp(-logit))
  d[ ,'Y'] <- rbinom(n, size = 1, prob = pY)

  return (d)
}

# TODO hyper params vs. params
setSigma_B_X_AH <- function(params, val) {
  params$sigma_B_X_AH <- val
  return(params)
}


setA <- function(params, val) {
  params$A <- val
  return(params)
}


set_sigma_B_X_AH <- function(params, val) {
  params$sigma_B_X_AH <- val
  return(params)
}
