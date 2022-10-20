# LearningWithExternalStats
Machine learning that combines internal data samples with limited statistics from external datasets.

## Installation

In addition to the main package we install `plpDataAdapter` to allow running all examples.
```
devtools::install_github("KI-Research-Institute/plpDataAdapter")
devtools::install_github("KI-Research-Institute/LearningWithExternalStats", build_vignettes = T)
```

## Basic example

A basic example with simulated data can be found using
```
vignette('logisticReweighting', package='LearningWithExternalStats')
```

## Example vignette for evaluation with OHDSI patient level prediction data

The following command opens a pdf vignette about using the package with 
[HADES PLP](https://ohdsi.github.io/PatientLevelPrediction/) data  
```
vignette('evaluatingPerformanceWithExternalStatsCDMExample', package='LearningWithExternalStats')
```
