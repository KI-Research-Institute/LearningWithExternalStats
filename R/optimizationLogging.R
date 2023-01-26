plotOptimizationLog <- function(l, wOptimizer, n) {

  batchSize <- min(wOptimizer$batchSize, n)
  batchesPerProbe <- wOptimizer$nProbe
  if (is.null(batchesPerProbe) || is.na(batchesPerProbe))
    batchesPerProbe <- 1
  nProbes <- nrow(l) / batchesPerProbe
  probeIdx <- (0:nProbes)*batchesPerProbe

  x <- 1:nrow(l)
  k <- ncol(l)
  par(mfrow=c(k,1))
  for (j in 1:k) {
    m <- colnames(l)[j]
    idx <- !is.na(l[,j])
    if (m == 'Primal objective')
      ylim <- c(0, mean(l[idx,j]+4*sd(l[idx,j])))
    else
      ylim <- NULL
    plot(x[idx], l[idx,j], type='l', xlab='Iteration', ylab = m, col='blue', ylim=ylim, cex=2, cex.lab=2, cex.axis=2)
    if (nProbes>1) {
      probeScore <- vector(mode = 'numeric', length = nProbes)
      for (probe in 1:nProbes)
        probeScore[probe] <- mean(l[(probeIdx[probe]+1):probeIdx[probe+1], j])
      lines(probeIdx[2:(nProbes+1)], probeScore, col='magenta', lw=2)
    }
  }
}
