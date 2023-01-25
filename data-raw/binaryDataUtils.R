normToBinary <- function(xy, freqs) {
  newxy <- data.frame(matrix(ncol = length(freqs)+1, nrow = nrow(xy)))
  colnames(newxy) <- c(colnames(xy)[1:length(freqs)], 'Y')
  for (i in 1:length(freqs)) {
    p <- pnorm(xy[, i], mean(xy[, i]) ,sd(xy[, i]))
    newxy[, i] <- as.numeric(p < freqs[i])
  }
  newxy[, 'Y'] <- xy[, 'Y']
  return(newxy)
}


normToBinaryTestData <- function(d, freqs=NULL) {
  xFeatures <- colnames(d$internalTest)[1:(ncol(d$internalTest)-1)]
  if (is.null(freqs)) {
    freqs = rbeta(length(xFeatures), 1, 10)
    cat(freqs, '\n')
  }

  d$internalTest <- normToBinary(d$internalTest, freqs)
  d$internalTrain <- normToBinary(d$internalTrain, freqs)
  d$externalTest <- normToBinary(d$externalTest, freqs)
  return(d)
}
