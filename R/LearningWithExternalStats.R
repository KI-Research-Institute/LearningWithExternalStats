# @file LearningWithExternalStats.R
#
# This file is part of LearningWithExternalStats.R
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' LearningWithExternalStats
#'
#' @description A package for running predictions using data in the OMOP CDM
#'
#' @docType package
#' @name LearningWithExternalStats
#' @import ParallelLogger
NULL

#' Internal train
#' @docType data
#' @keywords datasets
#' @name dIntTrain
#' @format A data frame containing the following elements:
#' \describe{
#'   \item{X1, ... , X10}{features}
#'   \item{Y}{outcome}
#' }
#' @usage
#' load(file=system.file('data/internalTrain.RData', package = "LearningWithExternalStats"))
#' TODO change this
NULL
