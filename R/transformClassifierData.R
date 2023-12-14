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
    d, transformType = 'Table 1', interactionVars = NULL, outcomeCol = 'Y', outcomeBalance = TRUE) {

  if (transformType == 'Table 1') {
    dTransformed <- computeTable1LikeTransformation(d, outcomeBalance=outcomeBalance, outcomeCol = outcomeCol)
  } else {
    rTransform <- switch(
      transformType,
      Interaction = reweightTransfrom$new(outcomeCol = outcomeCol, interactionVars = interactionVars),
      "Interaction v0" = reweightTransfrom$new(outcomeCol = outcomeCol, interactionVars = interactionVars, ver=0),
      reweightTransfrom$new(outcomeCol = outcomeCol)  # Default
    )
    cat('\n', transformType, 'formula:\n')
    f <- rTransform$getFormula(d)
    print(f)
    dTransformed <- model.matrix(f, data = d)
  }
  return(dTransformed)
}
