#' Post weighting diagnostics
#'
#' @description compute diagnostics of weighted samples
#'
#'
#' @param w a vector of weights
#' @param z a data frame of transformed feature-outcome pairs
#' @param mu a vector of means of transformed feature-outcome pairs
#'
#' @return a named list with the following:
#' \itemize{
#' \item {\code{Max weight}:} {Maximal weight}
#' \item {\code{chi2ToUnif}:} {Chi squared divergence between weights and a uniform distribution}
#' \item {\code{kl}:} {KL distance to uniform distribution}
#' \item {\code{Max Weighted SMD}:} {Maximum over features of SMD between internal set and external statistics}
#' }
#'
#' @export
postDiagnostics <- function(w, z, mu) {
  n <- length(w)
  p <- w/sum(w)
  klIdx <- p>0
  kl <- log(n) + as.numeric(t(p[klIdx]) %*% log(p[klIdx]))
  chi2ToUnif <- n*sum((p-1/n)**2)  #  = \sum(p-1/n)^2/1/n
  maxWeightedSMD <- computeMaxSMD(mu, z, p)
  ParallelLogger::logInfo(glue('Max WSMD {maxWeightedSMD}'))

  diagnostics <- list(
    'Max weight' = max(w),
    'chi2 to uniform' = chi2ToUnif,
    'kl' = kl,
    'Max Weighted SMD' = maxWeightedSMD)  # TODO add diagnostics
  return(diagnostics)
}

#' Get weighting results field names
#'
#' @return a character vector of field names
#'
#' @export
getWeightingResultsFieldNames <- function() {
  coreFields <- c('Opt err', 'n iter', 'n outcome',
    # 'Max weight', 'chi2 to uniform', 'kl',
    'Max Weighted SMD')
  fields <- coreFields
  for (f in fields) {
    fields <- c(fields, glue('95% lower {f}'), glue('Median    {f}'), glue('95% upper {f}'), glue('{f} sd'))
  }
  fields <- c('n repetitions', fields)
  return(fields)
}
