wfit <- function(trainer, X, Y, w = NULL) {
  if (missing(w))
    m <- trainer$fitFunction(X, Y, param = trainer$param)
  else
    m <- trainer$fitFunction(X, Y, trainer$param, w)
}


wpredict <- function(model, X) {
  return (model$predictFunction(model$model, X))
}


wimportant <- function(model) {
  return (model$importanceFunction(model))
}