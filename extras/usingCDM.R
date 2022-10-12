library(LearningWithExternalStats)
library(Eunomia)  # v1.0.2
library(FeatureExtraction)  # v3.2.0
library(PatientLevelPrediction)  # v6.0.4
library(Rmpfr); library(Matrix); library(dbplyr)  # To avoid Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'?

connectionDetails <- Eunomia::getEunomiaConnectionDetails()
Eunomia::createCohorts(connectionDetails)

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
  targetId = 1,
  outcomeDatabaseSchema = "main",
  outcomeTable = "cohort",
  outcomeIds = 3,
  cdmVersion = 5  #?
)

restrictPlpDataSettings <- createRestrictPlpDataSettings()  # sampleSize = 1000

plpData <- getPlpData(
    databaseDetails = databaseDetails,
    covariateSettings = covSettings,
    restrictPlpDataSettings = restrictPlpDataSettings)

summary(plpData)

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

lrResults <- runPlp(
    plpData = plpData,
    outcomeId = 3,
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

s <- lrResults$performanceEvaluation$evaluationStatistics
cat('Test AUROC:', s[(s['metric']=='AUROC') & (s['evaluation']=='Test'), 'value'][[1]], '\n')


cat('Data size =', length(lrResults$prediction$evaluationType),'\n')
for (s in unique(lrResults$prediction$evaluationType))
  cat(s, '\t', sum(lrResults$prediction$evaluationType==s), '\n')
tstIdx <- lrResults$prediction$evaluationType == 'Test'
cat('Outcome proportion:       ', mean(lrResults$prediction[tstIdx, 'outcomeCount']), '\n')
cat('Mean predicted proportion:', mean(lrResults$prediction[tstIdx, 'value']), '\n')  # predicted probability

externalDatabaseDetailes <- createDatabaseDetails(
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

resultsExternal <- externalValidateDbPlp(
  lrResults$model,
  validationDatabaseDetails = externalDatabaseDetailes)
