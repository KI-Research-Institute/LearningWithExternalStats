# @file WeightedEvaluationSummary.R
#

#' Compute weighted log-likelihood
#' 
#' @description compute weighted log-likelihood of a binary classifier predictions with respect to true labels 
#' 
#' @param Y a vector of binary labels
#' @param p a vector of numeric predictions (assuming in the range [0,1])
#' @param w a vector of non-negative weights
#' 
#' @return weighted log-likelihood
#' 
#' @export
WeightedLogLike <- function(Y, p, w=NULL) {
  if (is.factor(Y)) 
    Y <- as.numeric(Y)-1
  # handle zero likelihood instances
  idx0 <- (Y==0 & p==0) | (Y==1 & p==1)  # TODO consider a slack here
  llVec <- Y*log(p)+(1-Y)*log(1-p)
  llVec[idx0] <- 0  # using 0 log 0 = 0
  if (is.null(w))
    return(mean(llVec))
  else {
    return((llVec %*% w)/sum(w))
  }
}


WeightedBrier <- function(Y, p, w=NULL) {
  if (is.factor(Y)) 
    Yn <- as.numeric(Y)-1
  else
    Yn <- Y
  d <- (Yn-p)^2
  if (is.null(w))
    return(mean(d))
  else
    return((d %*% w)/sum(w))
}


meanPredictionRisk <- function(p, w=NULL) {
  if (is.null(w))
    return(mean(p))
  else
    return((p %*% w)/sum(w))
}


meanObservedRisk <- function(Y, w=NULL) {
  if (is.factor(Y)) 
    Yn <- as.numeric(Y)-1
  else
    Yn <- Y
  if (is.null(w))
    return(mean(Y))
  else
    return((Y %*% w)/sum(w))
}