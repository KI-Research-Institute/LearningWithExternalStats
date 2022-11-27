wrpart <- function() {
  param <- list()
  
  result <- list(
    name ='glmnet',
    param = param,
    fitFunction = wfitRpart
  )
  class(result) <- 'wmodel'
  return(result)
}


wfitRpart <- function(X, Y, param, w = NULL) {
  # todo change the API to formula?
  c <- rpart.control(maxdepth = 3, cp = 0)
  if (!(missing(w) || is.null(w))) {
    cat('Fitting weighted rpart\n')
    print(summary(w*length(w)))
  }
  d <- data.frame(cbind(X, Y))
  xFeatures <- glue('X{1:ncol(X)}')
  names(d) <- c(xFeatures, 'Y')
  formula <- paste(c('Y ~ 1', xFeatures), collapse = ' + ')
  if (missing(w))
    m <- rpart(as.formula(formula), data = d, method = 'class', control = c)
  else
    m <- rpart(as.formula(formula), data = d, weights = w, method = 'class', control = c)
  result <- list(model=m, predictFunction=wpredictRpart, importanceFunction=wimportantRpart, param=param)
  return(result)
}

wpredictRpart <- function(m, X) {
  d <- data.frame(X)
  xFeatures <- glue('X{1:ncol(X)}')
  names(d) <- xFeatures
  p <- predict(m, d)[,1]
  return(p)
}


wimportantRpart <- function(m) {
  return(NULL)
}


wprintRpart <- function(m) {
  print(m$model)
}
