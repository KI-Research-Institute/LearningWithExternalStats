# Generate a test model for unit tests and vignettes

rm(list=ls())
library(glmnet)
library(glue)

library(LearningWithExternalStats)

#  Load data
d <- LearningWithExternalStats::binaryAnchorData1

# Train and predict in internal and external sets
m <- ncol(d$internalTrain)-1
xTrain <- as.matrix(d$internalTrain[1:m])
yTrain <- as.vector(d$internalTrain[m+1])[[1]]
model1 <- cv.glmnet(xTrain, yTrain, family = "binomial", type.measure = "auc", alpha = 0.5)

usethis::use_data(binaryAnchorLR1, overwrite = T)
