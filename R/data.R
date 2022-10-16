#' Anchor data 1
#' @docType data
#' @keywords datasets
#' @name anchorData1
#' @format A named list containing the following elements:
#' \describe{
#'   \item{internalTrain}{A data frame of internal training set}
#'   \item{internalTest}{A data frame of internal test set}
#'   \item{externalTest}{A data frame of external test set}
#' }
#' Each dataset contains the following:
#' \describe{
#'   \item{X1, ... , X10}{features}
#'   \item{Y}{outcome}
#' }
#' @usage
#' anchorData1
NULL

#' Anchor logistic regression model 1
#' @docType data
#' @keywords datasets
#' @name anchorLR1
#' @format A logistic regression model trained using anchorData1::internalTrain
#' @usage
#' anchorLR1
NULL
