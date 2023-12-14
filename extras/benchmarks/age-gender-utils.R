library(glue)
library(reshape2)
library(R6)


#' Sample from table 1
#'
#' TODO simplify data.frame to vector
sampleAgeGender <- function(df, oColName, oValue, oSkew = NULL) {
  genderName <- 'gender = MALE'
  n <- df['n', oColName]
  ageRowsIndicator <- sapply(rownames(df), function(x) {substr(x, 0, 9) =='age group'})
  assertthat::assert_that(sum(ageRowsIndicator) > 0, msg = 'zero age rows in table')
  ageRows <- names(ageRowsIndicator)[ageRowsIndicator]
  agesFreqs <- df[ageRows, oColName]/100
  assertthat::assert_that(!any(is.na(agesFreqs)), msg = 'NA in ages frequencies')
  maleFreq <- df[genderName, oColName]/100

  clinicalRows <- rownames(df)[!(rownames(df) %in% c(genderName, 'n', '%', ageRows))]
  cat('n clinical features', length(clinicalRows), '\n')

  nAge <- length(ageRows)

  skewVec <- data.frame(rep(1, length(ageRows) + length(clinicalRows)))
  rownames(skewVec) <- c(ageRows, clinicalRows)

  if (!is.null(oSkew)) {
    inputSkewAgeFields <- intersect(ageRows, rownames(oSkew))
    inputSkewClinicalFields <- intersect(clinicalRows, rownames(oSkew))
    skewVec[inputSkewAgeFields, 1] <- oSkew[inputSkewAgeFields, 1]
    skewVec[inputSkewClinicalFields, 1] <- oSkew[inputSkewClinicalFields, 1]
  }
  cat('skews' ,skewVec[,1], '\n')

  # Sample gender column
  gender <- rbinom(n, 1, maleFreq)

  # Sample age columns given gender
  ageGroupIndicators <- matrix(nrow = n, ncol = length(agesFreqs))
  colnames(ageGroupIndicators) <- ageRows
  gSkew <- skewVec[ageRows, 1]
  probs0 <- agesFreqs / (1- maleFreq + gSkew * maleFreq)
  probs1 <- agesFreqs / (maleFreq + (1/gSkew) * (1 - maleFreq))
  nGender0 <- sum(gender == 0)
  nGender1 <- sum(gender == 1)
  if (nGender0 > 0)
    ageGroupIndicators[gender == 0, ] <- t(rmultinom(nGender0, 1, probs0))
  if (nGender1 > 1)
    ageGroupIndicators[gender == 1, ] <- t(rmultinom(nGender1, 1, probs1))

  Y <- rep(oValue, n)

  # Sample other columns if exist given gender
  if (length(clinicalRows) == 0) {
    d <- data.frame(cbind(gender, ageGroupIndicators[ ,1:(nAge-1)], Y))
  }
  else {
    clinicalFreqs <- df[clinicalRows, oColName]/100
    clinicalIndicators <- matrix(nrow = n, ncol = length(clinicalRows))
    colnames(clinicalIndicators) <- clinicalRows
    for (i in 1:length(clinicalRows)) {
      f <- clinicalRows[i]
      p0 <- clinicalFreqs[i] / (1- maleFreq + skewVec[f, 1] * maleFreq)
      p1 <- clinicalFreqs[i] / (maleFreq + (1/skewVec[f, 1]) * (1 - maleFreq))
      if (nGender0 > 0)
        clinicalIndicators[gender == 0, i] <- rbinom(nGender0, 1, p0)
      if (nGender1 > 0)
        clinicalIndicators[gender == 1, i] <- rbinom(nGender1, 1, p1)
    }
    d <- data.frame(cbind(gender, ageGroupIndicators[ ,1:(nAge-1)], clinicalIndicators, Y))
  }
  return(d)
}


#' Sample age gender outcome
#'
sampleAgeGenderOutcome <- function(df, dbName, oSkew = NULL) {
  d0 <- sampleAgeGender(df, glue('WithNoOutcome_CovariateMean_{dbName}'), oValue = 0, oSkew['WithNoOutcome'])
  d1 <- sampleAgeGender(df, glue('WithOutcome_CovariateMean_{dbName}'), oValue = 1, oSkew['WithOutcome'])

  d <- rbind(d0,d1)
  d <- d[sample(nrow(d)), ]
  return(d)
}


#' Test Single Age Gender Outcome
#'
#'
testSingleAgeGenderOutcome <- function(
    table1s, intDbName, extDbName, reweightConfigs,
    nMaxTrainSample = 50000, trainer = wglmnet(), skew = NULL, nRepetitions=10, maxCores=10, outputDir = getwd()) {

  dTrain <- sampleAgeGenderOutcome(table1s, intDbName)  # Assuming no skew in train set
  nTrainSample <- min(nMaxTrainSample, nrow(dTrain))
  cat('Downsampling from', nrow(dTrain), ' to ', nTrainSample, '\n')
  dTrain <- dTrain[sample(nrow(dTrain), nTrainSample), ]
  cat('p train Y', mean(dTrain[['Y']]), '\n')
  xFeatures <- colnames(dTrain)[1:(ncol(dTrain)-1)]
  internalX <- sapply(dTrain[xFeatures], as.numeric)
  m <- wfit(trainer, internalX, dTrain[['Y']])

  dTest <- sampleAgeGenderOutcome(table1s, intDbName)
  cat('p test Y', mean(dTest[['Y']]), '\n')
  dExt <- sampleAgeGenderOutcome(table1s, extDbName, skew)
  cat('p ext Y', mean(dExt[['Y']]), '\n')

  phatInt <- wpredict(m, as.matrix(dTest[xFeatures]))
  print(summary(phatInt))
  aucTest <- pROC::auc(as.factor(dTest[['Y']]), phatInt, direction="<", quiet=TRUE)
  cat('Internal AUROC', aucTest, '\n')

  phatExt <- wpredict(m, as.matrix(dExt[xFeatures]))
  cat('External predictions\n')
  print(summary(phatExt))
  aucExt <- pROC::auc(as.factor(dExt[['Y']]), phatExt, direction="<", quiet=TRUE)
  cat('External AUROC', aucExt, '\n')

  realResults <- list('Int AUROC' = aucTest, 'Ext AUROC' = aucExt)

  estimatedResults <- vector(mode = 'list', length = length(reweightConfigs))

  for (paramIdx in 1:length(reweightConfigs)) {

    cfg <- reweightConfigs[[paramIdx]]

    cat('--- Transform type', cfg[[1]], '--------------------------------------------------------------------------\n')
    dTransformedInt <- transformClassifierData(dTest, cfg[[1]], interactionVars = 'gender')
    dTransformedExt <- transformClassifierData(dExt, cfg[[1]], interactionVars = 'gender')
    muExt <- colMeans(dTransformedExt)

    internalData <- list(
      z = dTransformedInt,
      p = phatInt,
      y = dTest[['Y']]
    )

    externalEstimatorSettings <- createExternalEstimatorSettings(
      reweightAlgorithm = cfg[[2]],
      nRepetitions = nRepetitions,
      outputDir = outputDir,
      maxCores = maxCores
    )

    estimatedResults[[paramIdx]] <- estimateExternalPerformanceFromStatistics(
      internalData = internalData,
      externalStats = muExt,
      externalEstimatorSettings = externalEstimatorSettings
    )
    estAuroc <- estimatedResults[[paramIdx]]$estimation['AUROC', 'value']
    cat('Estimated AUROC', estAuroc, '\n')
  }
  return(list(realResults = realResults, estimatedResults = estimatedResults))
}



arrangeAgeGenderEmulationResults <- function(allResults, reweightConfigs) {
  # Similar to analyze2307
  keys <- c('analysis', 'internalDatabase', 'externalDatabase')
  summaryColNames <- c(
    keys, 'type', 'value.ext', 'outcomeCount', 'value.eval', 'opt.err', 'Max.Weighted.SMD', 'Estimation.Time')
  rSummary <- matrix(nrow=0, ncol=length(summaryColNames))
  colnames(rSummary) <- summaryColNames


  for (b in 1:length(allResults)) {
    experimentName <- names(allResults)[b]
    resultsName <- str_split(substring(experimentName, 1, nchar(experimentName)), '-')  # TODO learn
    internalDb <- resultsName[[1]][1]
    externalDb <- resultsName[[1]][3]


    if (!is.null(allResults[[b]]$realResults$`Int AUROC`)) {
      intAUC <- allResults[[b]]$realResults$`Int AUROC`
      if (!is.null(allResults[[b]]$realResults$`Ext AUROC`)) {
        extAUC <- allResults[[b]]$realResults$`Ext AUROC`
        l <- c(analysisIdx, internalDb, externalDb, 'internal', extAUC, NA, intAUC, NA, NA, NA)
        rSummary <- rbind(rSummary, l)

        for (paramIdx in 1:length(reweightConfigs)) {
          cfg <- reweightConfigs[[paramIdx]]
          wr <- allResults[[b]]$estimatedResults[[paramIdx]]$weightingResults
          estAuroc <- allResults[[b]]$estimatedResults[[paramIdx]]$estimation['AUROC', 'value']
          optErr <- wr[rownames(wr)=='Opt err', ]
          wsmd <- wr[rownames(wr)=='Max Weighted SMD', ]
          nO <- wr[rownames(wr)=='n outcome', ]
          estimationTime <- allResults[[b]]$estimatedResults[[paramIdx]]$estimationTime
          if (!is.null(estAuroc)) {
            l <- c(
              analysisIdx, internalDb, externalDb, cfg[[3]], extAUC, nO, estAuroc, optErr, wsmd, estimationTime)
            rSummary <- rbind(rSummary, l)
          }
        }
      }
    }
  }
  rSummary <- data.frame(rSummary)
  rSummary[['value.ext']] <- as.numeric(rSummary[['value.ext']])
  rSummary[['value.eval']] <- as.numeric(rSummary[['value.eval']])
  rSummary[['diff']] <- abs(rSummary[['value.eval']]-rSummary[['value.ext']])
  rSummary[['opt.err']] <- log10(as.numeric(rSummary[['opt.err']]))
  rSummary[['Max.Weighted.SMD']] <- as.numeric(rSummary[['Max.Weighted.SMD']])
  rSummary[['Estimation.Time']] <- as.numeric(rSummary[['Estimation.Time']])

  levels = c('internal')
  for (paramIdx in 1:length(reweightConfigs)) {
    cfg <- reweightConfigs[[paramIdx]]
    levels <- c(levels, cfg[[3]])
  }
  print(levels)
  rSummary$type <- factor(rSummary$type, levels = levels,ordered = TRUE)
  return(rSummary)
}

getSkewSpecificTag <- function(prefix, trainer, analysisIdx, nMaxTrainSample, nr, skewName) {
  tag <- glue("{prefix}{trainer$name}-a-{analysisIdx}-nTr-{nMaxTrainSample}-nr-{nr}-{skewName}")
}


plotTable1Results <- function(workDir, prefix, skewNames, nMaxTrainSample, nr, analysisIdx) {
  aSummary <- data.frame()
  for (skewName in skewNames) {
    if (nchar(skewName)>0) {
      skew <- read.csv(file.path(workDir, 'input', glue('{skewName}.csv')))
    } else {
      skew <- NULL
    }
    tag <- getSkewSpecificTag(prefix, trainer, analysisIdx, nMaxTrainSample, nr, skewName)
    rdsFileName <- file.path(workDir, 'output', glue('{tag}-results.RDS'))
    allResults <- readRDS(rdsFileName)
    cat('Read', rdsFileName, '\n')
    rSummary <- arrangeAgeGenderEmulationResults(allResults, reweightConfigs)
    rSummary[['skew']] <- skewName
    aSummary <- rbind(aSummary, rSummary)
  }
  aSummary[['has.eval']] = as.numeric(!is.na(aSummary[['value.eval']]))
  tag2 <- glue("{prefix}{trainer$name}-a-{analysisIdx}-nTrain-{nMaxTrainSample}-nr-{nr}-all-skews")
  p <- aSummary %>%
    ggplot(aes(x=type, y=diff, color = type)) + # group=type,
    geom_boxplot() +
    ylab('AUROC difference') +
    facet_wrap(~skew, nrow = 1)
  ggsave(file.path(workDir, 'output', glue('{tag2}.png')), height = 4, width = 4*length(skewNames))

  p <- aSummary %>%
    ggplot(aes(x=type, y=opt.err, color = type)) + # group=type,
    geom_boxplot() +
    ylab('log10(optimization error)') +
    facet_wrap(~skew, nrow = 1)
  ggsave(file.path(workDir, 'output', glue('{tag2}-opt.err.png')), height = 4, width = 4*length(skewNames))

  p <- aSummary %>%
    ggplot(aes(x=type, y=Max.Weighted.SMD, color = type)) + # group=type,
    geom_boxplot() +
    ylab('Max weighted SMD') +
    facet_wrap(~skew, nrow = 1)
  ggsave(file.path(workDir, 'output', glue('{tag2}-Max.Weighted.SMD.png')), height = 4, width = 4*length(skewNames))

  p <- aSummary %>%
    ggplot(aes(x=type, y=Estimation.Time, color = type)) + # group=type,
    geom_boxplot() +
    ylab('Estimation Time') +
    facet_wrap(~skew, nrow = 1)
  ggsave(file.path(workDir, 'output', glue('{tag2}-Estimation.Time.png')), height = 4, width = 4*length(skewNames))

  p <- aSummary %>%
    ggplot(aes(x=type, y=has.eval, color = type)) + # group=type,
    geom_boxplot() +
    ylab('Num evaluation') +
    facet_wrap(~skew, nrow = 1)
  ggsave(file.path(workDir, 'output', glue('{tag2}-has.eval.png')), height = 4, width = 4*length(skewNames))
}
