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
  rownames(ipdata) <- dataObject$labels$rowId
  return(ipdata)

}

makeFriendlyNames <- function(columnNames){

  columnNames <- gsub("[[:punct:]]", " ", columnNames)
  columnNames <- gsub(" ", "_", columnNames)
  return(columnNames)

}


summarizeResults <- function(s, evaluation) {

  cat(evaluation, 'AUROC:\t', s[(s['metric']=='AUROC') & (s['evaluation']==evaluation), 'value'][[1]], '\n')
  f <- c(
    'populationSize',
    'outcomeCount',
    'AUROC',
    '95% lower AUROC',
    '95% upper AUROC',
    # AUPRC
    'brier score',
    # brier score scaled
    # Eavg
    # E90
    # Emax
    'calibrationInLarge mean prediction',
    'calibrationInLarge observed risk'
  )
  print(format(s[(s['evaluation']==evaluation) & (s[['metric']] %in% f),c('metric', 'value')], digits=2))

}
