library(glue)


plotHighDimResults <- function(params) {
  testname <- getTestName(params)
  outputDir <- params$outputDir
  k <- length(params$estimationParams)
  legendLabels <- vector(mode = 'list', length = k)
  for (i in 1:k)
    legendLabels[[i]] <- params$estimationParams[[i]]$shortName
  plotRunTime(outputDir, testname, legendLabels)
  plotResults(outputDir, testname, legendLabels)
}


plotResults <- function(outputDir, testname, legendLabels) {
  k <- length(legendLabels)
  r = read.csv(file.path(outputDir, glue('{testname}.csv')))
  # print(r)
  colors <- c('red', 'purple', 'orange', 'blue', 'brown')
  pchs <- c(18,19,20,21,22)

  for (m in c('AUC', 'Brier', 'Calibration.prediction', 'Calibration.observed')) {
    allAUCs = c(r[[glue('Internal.{m}')]], r[[glue('External.{m}')]] ,
                r[[glue('Est..{m}.1')]] ,r[[glue('Est..{m}.2')]])
    idx <- !is.na(allAUCs)
    ylim = c(min(allAUCs[idx]), max(allAUCs[idx]))
    xlim = c(min(r[[glue('External.{m}')]]), max(r[[glue('External.{m}')]]))

    png(filename = file.path(outputDir, glue('{testname} {m}.png')))
    plot(r[[glue('External.{m}')]], r[[glue('Internal.{m}')]], ylim = ylim, pch=17,
         xlab = glue('External {m}'), ylab=glue('Estimation'))
    lwd = 1
    for (i in 1:k) {
      extResultsI <- r[[glue('External.{m}')]]
      points(extResultsI, r[[glue('Est..{m}.{i}')]], pch=pchs[i], col = colors[i])
      estAuc975 <- r[[glue('Est..{m}.{i}.high')]]
      estAuc025 <- r[[glue('Est..{m}.{i}.low')]]
      segments(extResultsI, estAuc025, extResultsI, estAuc975, col=colors[i], lwd = lwd, lty = 3);
    }
    lines(xlim, xlim)
    legend('bottomright', legend = c('Internal test', legendLabels), pch=c(17,pchs[1:k]),
           col=c('black', colors[1:k]))
    dev.off()
  }
}


plotRunTime <- function(outputDir, testname, legendLabels) {

  r = read.csv(file.path(outputDir, glue('{testname}.csv')))
  # print(r)
  colors <- c('red', 'purple', 'orange', 'blue', 'brown')
  n <- nrow(r)
  k <- length(legendLabels)
  runTimes <- matrix(ncol = k, nrow = n)

  for (i in 1:k) {
    c <- r[[glue('Estimation.Time.{i}')]]
    runTimes[, i] <- c
  }

  png(filename = file.path(outputDir, glue('{testname} runtimes.png')), width = 480, height = 480)
  boxplot(runTimes, col = colors[1:k], ylab = 'Runtime (minutes)')
  legend('bottomright', fill = colors[1:k], legend=legendLabels, cex=1, inset=0.02)
  dev.off()

}
