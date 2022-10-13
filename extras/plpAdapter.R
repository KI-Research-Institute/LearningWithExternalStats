transformPlpDataToDataFrame <- function(plpData, populationSettings, outcomeId) {
  #### new
  if(!is.null(plpData)){
    labels <- PatientLevelPrediction::createStudyPopulation(
      plpData = plpData,
      outcomeId = outcomeId,
      populationSettings = populationSettings
    )
  }
  # convert to matrix
  dataObject <- PatientLevelPrediction::toSparseM(
    plpData = plpData,
    cohort = labels
  )
  #sparse matrix: dataObject$dataMatrix
  #labels: dataObject$labels
  columnDetails <- as.data.frame(dataObject$covariateRef)

  cnames <- columnDetails$covariateName[order(columnDetails$columnId)]

  ipMat <- as.matrix(dataObject$dataMatrix)
  ipdata <- as.data.frame(ipMat)
  colnames(ipdata) <-  makeFriendlyNames(cnames)
  ipdata$outcome <- dataObject$labels$outcomeCount
  return(ipdata)
}

makeFriendlyNames <- function(columnNames){

  columnNames <- gsub("[[:punct:]]", " ", columnNames)
  columnNames <- gsub(" ", "_", columnNames)
  return(columnNames)

}
