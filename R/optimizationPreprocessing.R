#' Normalize data and expectations
#'
#' @description
#'
#' Normalize data and expectations.
#'
#' @param Z data frame
#' @param mu vector of expectations
#' @param minSd float value of minimum standard deviation. Columns with a smaller sd will be removed
#'
#' @return
#' list
#'
#' @export
normalizeDataAndExpectations <- function(Z, mu, minSd) {
  # Preprocess and filter columns with low variance
  muZ <- colMeans(Z)
  sdZ <- apply(Z, 2, sd)
  useCols <- sdZ>minSd # remove columns with low variance
  if (sum(!useCols) > 0) {
    ParallelLogger::logInfo(
      sprintf("Removing %d columns with low variance (<%.2g). Please examing this columns.", sum(!useCols), minSd))
    Z <- Z[,names(Z)[useCols]]  # remove low variance columns todo: flag a positivity issue
    sdZ <- sdZ[useCols]
    muZ <- muZ[useCols]
    mu <- mu[useCols]
  }
  # Normalize
  p <- ncol(Z)
  for (i in 1:p) Z[,i] <- (Z[,i]-muZ[i])/sdZ[i]
  mu <- (mu-muZ)/sdZ
  return(list(Z=Z, mu=mu, muZ=muZ, sdZ=sdZ))
}
