rm(list=ls())
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./anchorModel.R')

n <- 1000
hyperParams = getDefaultAnchorHyperParams()
hyperParams$sigma_B_X_AH = 0.5  # Medium shift
params <- anchorParams(hyperParams)
params$A <- 0
internalTrain <- sampleAnchor(params, n)
internalTest <- sampleAnchor(params, n)
params$A <- 1
externalTest <- sampleAnchor(params, n)

anchorData1 <- list(internalTrain = internalTrain, internalTest = internalTest, externalTest = externalTest)

usethis::use_data(anchorData1, overwrite = TRUE)
# Run install()
