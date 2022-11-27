library(glue)

plotHighDimResults <- function(params) {
  testname <- getTestName(params)
  outputDir <- params$outputDir
  r = read.csv(file.path(outputDir, glue('{testname}.csv')))
  # print(r)
  allAUCs = c(r[['Internal.AUC']], r[['External.AUC']] ,r[['Est..Ext..AUC.1']] ,r[['Est..Ext..AUC.2']])
  idx <- !is.na(allAUCs)
  ylim = c(min(allAUCs[idx]), max(allAUCs[idx]))
  xlim = c(min(r[['External.AUC']]), max(r[['External.AUC']]))

  png(filename = file.path(outputDir, glue('{testname}.png')))
  plot(r[['External.AUC']], r[['Internal.AUC']], ylim = ylim, pch=17, xlab = 'External AUC', ylab='Estimation')
  points(r[['External.AUC']], r[['Est..Ext..AUC.1']], pch=18, col = 'red')
  points(r[['External.AUC']], r[['Est..Ext..AUC.2']], pch=19, col = 'purple')
  lines(xlim, xlim)
  legend('bottomright', c('Internal test', 'Estimation 1', 'Estimation 2'), pch=c(17,18,19),
         col=c('black', 'red', 'purple'))
  dev.off()
}
