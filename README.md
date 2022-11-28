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

### Pre-weighting diagnostics
This stage tests the feasibility of re-weighting the internal set to achieve external statistics. 

### Re-weighting
This stage re-weights the internal data such that the weighted set has similar statistics to the external ones.
The internal re-weighting algorithm can handle up to a few 10,000's samples and up to 500 features.
In case the number of samples is large, the users can determine a maximum number of samples and the evaluation
algorithm will invoke internal re-weighting using sub-samples of this size (in this case we recommend 5000-10,000).

### Post-weighting diagnostics
If the re-weighting algorithm finds a vector of weights that satisfies the task's constraints, the evaluation module
further examines them. The main examination is computation of the maximum weighted *standardized mean difference* (SMD).
The algorithm continues if the maximum weighted SMD is greater than 0.2.

### Estimation
Finally, performance measures are computed using the weighted versions.
