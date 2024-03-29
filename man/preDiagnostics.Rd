% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preDiagnostics.R
\name{preDiagnostics}
\alias{preDiagnostics}
\title{Pre re-weighting diagnostics}
\usage{
preDiagnostics(z, mu, maxDiff, npMinRatio = 4, maxSubset = 20000)
}
\arguments{
\item{z}{a data frame of transformed feature-outcome pairs}

\item{mu}{a vector of means of transformed feature-outcome pairs}

\item{maxDiff}{maximum difference from which to check max prop}

\item{npMinRatio}{minimum ratio between number of features and number of examples}

\item{maxSubset=20000}{}
}
\value{
a named list with the following fields:
\itemize{
\item{outOfRange}{list of external statistics is out of range of a numeric feature}
\item{representedFeatures}{names of features }
\item{zidx}{indices of samples that should be used in reweighting}
\item{highlySkewedBinary}{list of highly skewed variables}
\item{status}{a string that determines}
\item{structuredLog}{A data frame with indicators of prediagnosis results and descriptions, see below}
}

Structured log is a data-frame with a single value column and the following row names:
\describe{
\item{preDiagnosisStatus}{over all status of diagnostics}
\item{missingInMu}{there are variables in the internal data matrix that are missing in the external statistics}
\item{missingInZ}{variable in the external statistics are missing in the internal data}
\item{nOverlapFeatures}{number of overlapping features}
\item{naInZ}{NA values in the internal data}
\item{naInMu}{NA values in the statistics vector}
\item{highlySkewedBinary}{highly skewed binary variables, generally indicating a few sample in the internal
dataset but a high proportion of them in the external ones}
\item{highSkewRepresentative'}{an example for a highly skewed variable}
\item{incompatableUnaryVariable'}{value of an internal unary variable is different than the external one more than
a threshold}
\item{incompatibleUnaryRepresentative'}{example for incompatible unary variable}
\item{outOfRange}{externa statistic is out of range}
\item{outOfRangeRepresentative}{example}
\item{fewSamples}{few sample relative to number of }
\item{fewSamplesDescription}{elaboration on number of samples and dimension}
}
}
\description{
compare external expectations and internal means before reweighting

The following tests are made:
\enumerate{
\item z and mu have the same set of variables
\item z and mu do not contain NA entries
\item skew in binary features
}

Additionally, internal samples are removed when ...
}
