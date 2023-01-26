wglmnet <- function(alpha=0.5, ntop=1000) {
  param <- list(alpha=alpha, ntop=ntop)

  result <- list(
    name ='glmnet',
    param = param,
    fitFunction = wfitGlmnet
  )
  class(result) <- 'wmodel'
  return(result)
}


wfitGlmnet <- function(X, Y, param, w = NULL, verbose = F) {
  # todo change the API to formula?
  if (verbose)
    if (!(missing(w) || is.null(w))) {
      cat('Fitting weighted Glmnet\n')
      print(summary(w*length(w)))
    }
  if (missing(w))
    m <- cv.glmnet(X, Y, family = "binomial", type.measure = "auc", alpha = param$alpha)
  else
    m <- cv.glmnet(X, Y, weights = w, family = "binomial", type.measure = "auc", alpha = param$alpha)
  result <- list(model=m, predictFunction=wpredictGlmnet, importanceFunction=wimportantGlmnet, param=param)
  return(result)
}

wpredictGlmnet <- function(m, X) {
  p <- predict(m, X, type = "response", s = "lambda.min")[,1]  # 1se
  return(p)
}


wimportantGlmnet <- function(m) {
  # TODO - add additional info if exists
  beta <- coef(m$model, s = "lambda.min")
  betaDf <- data.frame(name = beta@Dimnames[[1]][beta@i + 1], coefficient = beta@x)
  nr <- nrow(betaDf)
  ntop <- min(m$param$ntop, nr-1)
  betaDf <- betaDf[2:nr, ]
  betaDf <- betaDf[order(-abs(betaDf$coefficient))[1:ntop], ]
  return(betaDf[, 'name'])
}


wprintGlmnet <- function(m) {
  beta <- coef(m$model, s = "lambda.min")
  betaDf <- data.frame(name = beta@Dimnames[[1]][beta@i + 1], coefficient = beta@x)
  cat('GLMnet coefficients:\n')
  print(betaDf)
}
