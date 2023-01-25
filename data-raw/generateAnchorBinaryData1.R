rm(list=ls())
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
source('./binaryDataUtils.R')

d <- LearningWithExternalStats::anchorData1
binaryAnchorData1 <- normToBinaryTestData(d)

usethis::use_data(binaryAnchorData1, overwrite = TRUE)

