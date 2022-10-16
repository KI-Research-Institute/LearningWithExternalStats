# Generate a test model for unit tests and vignettes

rm(list=ls())
library(glmnet)
library(glue)

library(LearningWithExternalStats)

#  Load data
d <- LearningWithExternalStats::anchorData1

# Train and predict in internal and external sets
xFeatures <- colnames(d$internalTest)[1:(ncol(d$internalTest)-1)]
anchorLR1 <- cv.glmnet(sapply(d$internalTrain[xFeatures], as.numeric), d$internalTrain[['Y']],
               family = "binomial", type.measure = "auc", alpha = 0.5)

usethis::use_data(anchorLR1, overwrite = T)
