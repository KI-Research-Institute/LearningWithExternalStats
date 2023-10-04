transformByType <- function(transformType, d, interactionVar = NULL) {

  if (transformType == 'Current flat') {
    cat('Current flat\n')
    dTransformed <- computeTable1LikeTransformation(d, outcomeBalance=TRUE)
  } else {
    if (transformType == 'New flat') {
      rTransform <- reweightTransfrom$new(outcomeCol = 'Y')
    } else {
      rTransform <- reweightTransfrom$new(outcomeCol = 'Y', interactionVars = interactionVar)
    }
    cat('\n', transformType, 'formula:\n')
    f <- rTransform$getFormula(d)
    print(f)
    dTransformed <- model.matrix(f, data = d)
  }
  return(dTransformed)
}
