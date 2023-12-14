#' @importFrom R6 R6Class
#' @import glue
NULL


# @file reweightTransfroms.R
#
# Copyright 2021 KI Research Institute
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


#' @title a re-weight transform that supports interactions among features
#'
#' @description a class that generates transform formulas
#'
#' @field outcomeCol name of outcome column
#' @field interactionVars a character or character vector with second order interaction variables. For every variable
#' in this field, an interaction is added with all other variables.
#'
#' @export
reweightTransfrom <- R6::R6Class(
  "reweightTransfrom",

  public = list(

    # TODO consider adding a flag of outcomeBalance, conseider deprecating old interaction formula version
    outcomeCol = NA,
    interactionVars = NA,
    ver = 1,

    #' @param outcomeCol outcome column name
    #' @param interactionVars names of interaction variables (NULL, character or character vector)
    #' @param ver version of interaction formula. The most updated is 1
    #'
    initialize = function(outcomeCol = 'Y', interactionVars = NULL, ver=1) {
      self$outcomeCol <- outcomeCol
      self$interactionVars <- interactionVars
      self$ver <- ver
    },

    #' @param X a data frame in numeric format
    getFormula = function(X) {

      xTerms <- glue('(. - {self$outcomeCol})')
      # Initialize the formula with flat terms, in case there are variables with >2 values we add their squares
      gt2Values <- sapply(X, function(x) length(unique(x)) > 2)
      gt2Values <- gt2Values & (colnames(X) != self$outcomeCol)  # Remove outcome from squared features
      if (sum(gt2Values) > 0) {
        quadTerms <- paste(glue('I({colnames(X)[gt2Values]}^2)'), collapse = ' + ')
        singleVarTerms <- glue('(. + {quadTerms} - {self$outcomeCol})')
      }
      else
        singleVarTerms <- xTerms
      Table1Terms <- glue('I({self$outcomeCol}):{singleVarTerms}+I(1-{self$outcomeCol}):{singleVarTerms}')
      fs <- glue('~ -1 + {self$outcomeCol} + {Table1Terms}')

      # Add interaction if given
      if (!is.null(self$interactionVars)) {
        if (self$ver == 1) {
          xInteraction <- glue('(. - {self$outcomeCol} - {self$interactionVars[1]}):{self$interactionVars[1]}')
          if (length(self$interactionVars)>1) {
            for (i in 2:length(interactionVars)) {
              xi <- glue('(. - {self$outcomeCol} - {self$interactionVars[i]}):{self$interactionVars[i]}')
              xInteraction <- glue('{xInteraction} + {xi}')
            }
            xInteraction <- glue('({xInteration})')
          }
          fs <- glue('{fs} + I({self$outcomeCol}):{xInteraction}+I(1-{self$outcomeCol}):{xInteraction}')
        } else {
          interactionVarsTerm <- paste(self$interactionVars, collapse = ' + ')
          if (length(self$interactionVars) > 1)
            interactionVarsTerm <- glue("({interactionVarsTerm})")
          fs <- glue('{fs} + ({Table1Terms}):{interactionVarsTerm}')
        }
      }
      f <- as.formula(fs)
      return(f)

    }

  )
)





#' @title Compute Table 1 like transformation
#'
#' @description
#'
#' Compute features and squared numeric feature multiplied by outcome and 1-outcome.
#'
#' @param X data frame
#' @param outcomeBalance boolean specifying if to use outcome in balancing
#' @param outcomeCol string specifying the name of the outcome column
#'
#' @return
#' A data frame
#'
#' @export
computeTable1LikeTransformation <- function(X, outcomeBalance, outcomeCol='Y') {
  # Add squares of numeric features
  is_numeric <- sapply(X, function(x) length(unique(x))>2)  # TODO - consider a more elegant way
  for (s in names(is_numeric[is_numeric])) {
    sNew <- paste('(' ,s, '^2)', sep='')
    X[[sNew]] <- (X[[s]]**2)
  }
  # Convert binary CVXR::Variables to numeric 0-1
  is_factor <- sapply(X, is.factor)
  X[is_factor] <- sapply(X[is_factor], as.numeric)
  X[is_factor] <- X[is_factor] - 1
  # Add interactions with outcome
  if (outcomeBalance) {
    Z = data.frame(row.names = row.names(X))
    Z[['Y']] <- X[[outcomeCol]]  # TODO decide if to maintain this
    x_names <- names(X)
    ext_x <- x_names[x_names != outcomeCol]
    ext_x_y1 <- paste(ext_x, ':y1', sep="")
    for (i in 1:length(ext_x)) {
      Z[[ext_x_y1[i]]] <- X[[outcomeCol]] * X[[ext_x[i]]]
    }
    ext_x_y0 <- paste(ext_x, ':y0', sep="")
    for (i in 1:length(ext_x)) {
      Z[[ext_x_y0[i]]] <- (1-X[[outcomeCol]]) * X[[ext_x[i]]]
    }
  } else {
    Z <- X
  }
  return (Z)
}
