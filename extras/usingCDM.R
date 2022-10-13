rm(list = ls())

library(LearningWithExternalStats)
library(Eunomia)  # v1.0.2
library(FeatureExtraction)  # v3.2.0
library(PatientLevelPrediction)  # v6.0.4
library(Rmpfr); library(Matrix); library(dbplyr)  # To avoid Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'?
library(glue)

# TODO extract test samples
# TODO revive MaxSMD

source('./extras/plpAdapter.R')

# Definitions
outcomeId <- 3
internalTargetId <- 1  #
externalTargetId <- 2



summarizeResults <- function(s, evaluation) {
  cat(evaluation, 'AUROC:', s[(s['metric']=='AUROC') & (s['evaluation']==evaluation), 'value'][[1]], '\n')

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

  # cat('Data size =', length(lrResults$prediction$evaluationType),'\n')
  # for (s in unique(lrResults$prediction$evaluationType))
  #   cat(s, '\t', sum(lrResults$prediction$evaluationType==s), '\n')
  # cat('Outcome proportion:       ', mean(lrResults$prediction[tstIdx, 'outcomeCount']), '\n')
  # cat('Mean predicted proportion:', mean(lrResults$prediction[tstIdx, 'value']), '\n')  # predicted probability
}

# Activate Eunomia demo database
connectionDetails <- Eunomia::getEunomiaConnectionDetails()
Eunomia::createCohorts(connectionDetails)

# Study setup
covSettings <- createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useConditionGroupEraLongTerm = TRUE,
  useConditionGroupEraAnyTimePrior = TRUE,
  useDrugGroupEraLongTerm = TRUE,
  useDrugGroupEraAnyTimePrior = TRUE,
  useVisitConceptCountLongTerm = TRUE,
  longTermStartDays = -365,
  endDays = -1)

databaseDetails <- createDatabaseDetails(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "main",
  cdmDatabaseName = '',
  cdmDatabaseId = '',  # ?
  tempEmulationSchema = NULL,  # is this important to avoid further errors in getPlpData?
  cohortDatabaseSchema = "main",
  cohortTable = "cohort",
  targetId = internalTargetId,
  outcomeDatabaseSchema = "main",
  outcomeTable = "cohort",
  outcomeIds = outcomeId,
  cdmVersion = 5  #?
)

restrictPlpDataSettings <- createRestrictPlpDataSettings()  # sampleSize = 1000

# Generate data
plpData <- getPlpData(
    databaseDetails = databaseDetails,
    covariateSettings = covSettings,
    restrictPlpDataSettings = restrictPlpDataSettings)

summary(plpData)

# Prediction model definitions
populationSettings <- createStudyPopulationSettings(
  washoutPeriod = 364,
  firstExposureOnly = FALSE,
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback = 9999,
  riskWindowStart = 1,
  riskWindowEnd = 365,
  minTimeAtRisk = 364,
  requireTimeAtRisk = TRUE,
  includeAllOutcomes = TRUE)

splitSettings <- createDefaultSplitSetting(
  trainFraction = 0.75,
  testFraction = 0.25,
  type = 'stratified',
  nfold = 2,
  splitSeed = 1234
)

sampleSettings <- createSampleSettings()
featureEngineeringSettings <- createFeatureEngineeringSettings()
preprocessSettings <- createPreprocessSettings(
  minFraction = 0.01,
  normalize = T,
  removeRedundancy = T
)


lrModel <- setLassoLogisticRegression()

executeSettings = createExecuteSettings(
  runSplitData = T,
  runSampleData = T,
  runfeatureEngineering = T,
  runPreprocessData = T,
  runModelDevelopment = T,
  runCovariateSummary = T
)

# Train the model
lrResults <- runPlp(
    plpData = plpData,
    outcomeId = outcomeId,
    analysisId = 'singleDemo',
    analysisName = 'Demonstration of runPlp for training single PLP models',
    populationSettings = populationSettings,
    splitSettings = splitSettings,
    sampleSettings = sampleSettings,
    featureEngineeringSettings = featureEngineeringSettings,
    preprocessSettings = preprocessSettings,
    modelSettings = lrModel,
    logSettings = createLogSettings(),
    executeSettings = executeSettings,
    saveDirectory = file.path(getwd(), 'singlePlp')
  )



tstIdx <- lrResults$prediction$evaluationType == 'Test'
cvIdx <- lrResults$prediction$evaluationType == 'CV'
trainIdx <- lrResults$prediction$evaluationType == 'Train'

summarizeResults(lrResults$performanceEvaluation$evaluationStatistics, 'Test')

externalDatabaseDetailes <- createDatabaseDetails(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "main",
  cdmDatabaseName = '',
  cdmDatabaseId = '',  # ?
  tempEmulationSchema = NULL,  # is this important to avoid further errors in getPlpData?
  cohortDatabaseSchema = "main",
  cohortTable = "cohort",
  targetId = externalTargetId,
  outcomeDatabaseSchema = "main",
  outcomeTable = "cohort",
  outcomeIds = outcomeId,
  cdmVersion = 5  #?
)

resultsExternal <- externalValidateDbPlp(
  lrResults$model,
  validationDatabaseDetails = externalDatabaseDetailes)


trainData <- transformPlpDataToDataFrame(plpData, populationSettings, outcomeId = 3)


externalDatabaseDetails <- createDatabaseDetails(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "main",
  cdmDatabaseName = '',
  cdmDatabaseId = '',  # ?
  tempEmulationSchema = NULL,  # is this important to avoid further errors in getPlpData?
  cohortDatabaseSchema = "main",
  cohortTable = "cohort",
  targetId = 2,
  outcomeDatabaseSchema = "main",
  outcomeTable = "cohort",
  outcomeIds = 3,
  cdmVersion = 5  #?
)

externalPlpData <- getPlpData(
  databaseDetails = externalDatabaseDetails,
  covariateSettings = covSettings,
  restrictPlpDataSettings = restrictPlpDataSettings)

externalData <- transformPlpDataToDataFrame(externalPlpData, populationSettings, outcomeId = outcomeId)


covariateImportance <- lrResults$model$covariateImportance
importantCovariates <- covariateImportance[abs(covariateImportance['covariateValue']) > 0, ]
estimationCovariates <- importantCovariates[['covariateName']]
estimationCovariates <- makeFriendlyNames(estimationCovariates)

dim(trainData[estimationCovariates])
dim(externalData[estimationCovariates])

ZInt <- computeTable1LikeTransformation(
  trainData[c(estimationCovariates, 'outcome')], outcomeBalance = T, outcomeCol = 'outcome')
zext <- computeTable1LikeTransformation(
  externalData[c(estimationCovariates, 'outcome')], outcomeBalance = T, outcomeCol = 'outcome')

muExt <- colMeans(zext)

# TODO - check if we need to hand the outcome
# TODO add a test to the package to make sure no: Error in rep(TRUE, n1) : invalid 'times' argument. check internalData
# TODO check vector lengths

divergence <- 'entropy'
lambda <- 1e-2
minW <- 1e-6
optimizationMethod <- 'dual'

evalIdx <- tstIdx | trainIdx
prediction <- lrResults$prediction[evalIdx, ]
prediction <- prediction[order(prediction$rowId), ]

if (sum(prediction$outcomeCount != trainData[['outcome']]) > 0)  # TODO make if more rigorous using rowId
  stop('need to allign prediction and trainData')

internalData <- list(z=ZInt, p = prediction$value, y = trainData[['outcome']])
estimatedResults <- estimateExternalPerformanceFromStats(
  internalData, muExt, divergence = divergence, lambda = lambda, minW = minW, optimizationMethod = optimizationMethod,
  nboot = 0)
summarizeResults(lrResults$performanceEvaluation$evaluationStatistics, 'Test')
summarizeResults(resultsExternal[[2]]$performanceEvaluation$evaluationStatistics, 'Validation')

# TODO change into a data-frame
cat(glue('KL divergence between estimated weights and uniform ones = {format(estimatedResults$summary$kl, digits=3)}'),'\n')
cat(glue('Estimated AUC   = {format(estimatedResults$summary$AUC, digits=3)}'),'\n')
cat(glue('Estimated Brier = {format(estimatedResults$summary$Brier, digits=3)}'),'\n')
cat(glue('Predicted Risk  = {format(estimatedResults$summary$predictedR, digits=3)}'),'\n')
cat(glue('Estimated Risk  = {format(estimatedResults$summary$observedR, digits=3)}'),'\n')


