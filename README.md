# LearningWithExternalStats
Machine learning that combines internal data samples with limited statistics from external datasets.

## Installation

In addition to the main package we install `plpDataAdapter` to allow running examples
that use [HADES patient level prediction](https://ohdsi.github.io/PatientLevelPrediction/) data:
```
devtools::install_github("KI-Research-Institute/plpDataAdapter")
devtools::install_github("KI-Research-Institute/LearningWithExternalStats")
```


## Overview

This package estimates the performance of a predictive model of an external data-set when there is only access to
limited statistics in the external set and to detailed samples in the internal one. The main function works in
three stages:

### Pre-diagnostics
This stage tests the feasibility of re-weighting the internal set to achieve external statistics. 

### Re-weighting
This stage re-weights the internal data such that the weighted set has similar statistics to the external ones.

### Post-diagnostics


### Estimation
