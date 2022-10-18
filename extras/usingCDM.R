rm(list = ls())

library(LearningWithExternalStats)
library(Eunomia)  # v1.0.2
library(FeatureExtraction)  # v3.2.0
library(PatientLevelPrediction)  # v6.0.4
library(Rmpfr); library(Matrix); library(dbplyr)
# To avoid Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'?
library(glue)
source('./extras/plpAdapter.R')
# TODO extract test samples, Revive MaxSMD

# Definitions
outcomeId <- 3
internalTargetId <- 1
externalTargetId <- 2

# _____________________________________________________________________________________________________________________
#
# Activate Eunomia demo database
#
connectionDetails <- Eunomia::getEunomiaConnectionDetails()
Eunomia::createCohorts(connectionDetails)

# _____________________________________________________________________________________________________________________
#
# Study setup
#
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

# _____________________________________________________________________________________________________________________
#
# Internal node:
#

# Generate internal data
#
internalPlpData <- getPlpData(
    databaseDetails = databaseDetails,
    covariateSettings = covSettings,
    restrictPlpDataSettings = restrictPlpDataSettings)
summary(internalPlpData)

# Prediction model definitions
#
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
#
internalResults <- runPlp(
    plpData = internalPlpData,
    outcomeId = outcomeId,
    analysisId = '../singleDemo',  # TODO
    analysisName = 'Demonstration of runPlp for training single PLP models',
    populationSettings = populationSettings,
    splitSettings = splitSettings,
    sampleSettings = sampleSettings,
    featureEngineeringSettings = featureEngineeringSettings,
    preprocessSettings = preprocessSettings,
    modelSettings = lrModel,
    logSettings = createLogSettings(),
    executeSettings = executeSettings,
    saveDirectory = file.path(getwd(), '..', 'singlePlp') # TODO
  )

summarizeResults(internalResults$performanceEvaluation$evaluationStatistics, 'Test')

# Identify model covariates. These will be used in the estimation algorithm.
#
covariateImportance <- internalResults$model$covariateImportance
importantCovariates <- covariateImportance[abs(covariateImportance['covariateValue']) > 0, ]
estimationCovariates <- importantCovariates[['covariateName']]
estimationCovariates <- makeFriendlyNames(estimationCovariates)


# _____________________________________________________________________________________________________________________
#
# External node: generate data, evaluate model, transform and compute means
#
externalDatabaseDetails <- createDatabaseDetails(
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

externalPlpData <- getPlpData(
  databaseDetails = externalDatabaseDetails,
  covariateSettings = covSettings,
  restrictPlpDataSettings = restrictPlpDataSettings)

externalResults <- externalValidateDbPlp(  # TODO how should we avoid duplicated generation of plpData?
  internalResults$model,
  validationDatabaseDetails = externalDatabaseDetails)

externalXY <- transformPlpDataToDataFrame(externalPlpData, populationSettings, outcomeId = outcomeId)
# TODO verify that estimation covariates are included
ZExt <- computeTable1LikeTransformation(
  externalXY[c(estimationCovariates, 'outcome')], outcomeBalance = T, outcomeCol = 'outcome')
muExt <- colMeans(ZExt)

# _____________________________________________________________________________________________________________________
#
# Internal node: evaluation of external dataset peformence using internal dataset
#

# Identify and align predictions with internal data
tstIdx <- internalResults$prediction$evaluationType == 'Test'
cvIdx <- internalResults$prediction$evaluationType == 'CV'
trainIdx <- internalResults$prediction$evaluationType == 'Train'
evalIdx <- tstIdx | trainIdx
prediction <- internalResults$prediction[evalIdx, ]
prediction <- prediction[order(prediction$rowId), ]

# Prepare input for evaluation
internalXY <- transformPlpDataToDataFrame(internalPlpData, populationSettings, outcomeId = 3)
if (sum(prediction$outcomeCount != internalXY[['outcome']]) > 0)  # TODO make if more rigorous using rowId
  stop('need to allign prediction and internalXY')
ZInt <- computeTable1LikeTransformation(
  internalXY[c(estimationCovariates, 'outcome')], outcomeBalance = T, outcomeCol = 'outcome')

# TODO - check if we need to hand the outcome
# TODO add a test to the package to make sure no: Error in rep(TRUE, n1) : invalid 'times' argument. check internalData
# TODO check vector lengths

internalData <- list(z=ZInt, p = prediction$value, y = internalXY[['outcome']])
estimatedResults <- estimateExternalPerformanceFromStats(  # TODO automatic selection of parameters?
  internalData, muExt, divergence = 'entropy', lambda = 1e-2, minW = 1e-6, optimizationMethod = 'dual', nboot = 10)

# _____________________________________________________________________________________________________________________
#
# Compare evaluation and real results
#

summarizeResults(internalResults$performanceEvaluation$evaluationStatistics, 'Test')
summarizeResults(externalResults[[2]]$performanceEvaluation$evaluationStatistics, 'Validation')

cat('Estimated metrics:\n')
# TODO change into a data-frame
print(estimationSummaryToDF(estimatedResults$summary))
