#' @title Transform classifier data
#'
#' @description
#'
#' Transform features and outcome columns, for example concatenate squared features, feature-outcome interaction
#' and more.
#'
#' @param d data frame
#' @param transformType `'Table 1'` or `'Interaction'` or `'Flat'`
#' @param interactionVars interaction variables in case `transformType` is `'Interaction'`
#' @param outcomeBalance boolean specifying if to use outcome in balancing
#' @param outcomeCol string specifying the name of the outcome column
#'
#' @return
#' A data frame
#'
#' @export
transformClassifierData <- function(
    d, transformType = 'Table 1', interactionVar = NULL, outcomeCol = 'Y', outcomeBalance = outcomeBalance) {

  if (transformType == 'Table 1') {
    dTransformed <- computeTable1LikeTransformation(d, outcomeBalance=TRUE, outcomeCol = outcomeCol)
  } else {
    if (transformType == 'Interaction') {
      rTransform <- reweightTransfrom$new(outcomeCol = outcomeCol, interactionVars = interactionVar)
    } else {
      rTransform <- reweightTransfrom$new(outcomeCol = outcomeCol)
    }
    cat('\n', transformType, 'formula:\n')
    f <- rTransform$getFormula(d)
    print(f)
    dTransformed <- model.matrix(f, data = d)
  }
  return(dTransformed)
}
