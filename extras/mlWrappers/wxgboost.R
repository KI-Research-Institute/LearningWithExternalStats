require(xgboost)
library(mlr3)
library(mlr3pipelines)
library(mlr3verse)
library(mlr3tuning)
library(glue)
# TODO - generalize to mlr3

wXGBoost <- function(gamma=10, eta = 0.1, max_depth = p_int(lower = 2, upper = 6), nrounds = 50) {
  lgr::get_logger("mlr3")$set_threshold("warn")
  lgr::get_logger("bbotk")$set_threshold("warn")
  
  xgboostlearner <- lrn(
    "classif.xgboost", objective="binary:logistic", 
    eval_metric='logloss',  # 'error', 'logloss' 'auc'
    gamma = gamma, booster = "gbtree",
    eta = 0.1, max_depth = 6L, nrounds = 50,
    predict_type = "prob")
  # eta = eta, max_depth = max_depth, nrounds = nrounds,  # TODO - can we remove this
  
  search_space = ps(
    eta = p_dbl(lower = 0.001, upper = eta),
    max_depth = max_depth,
    nrounds = p_int(lower = 1, upper = nrounds)  # lower = 20
  )

  param <- list(
    methodType='mlr3',
    name='mlr3 classif.xgboost',
    learner=xgboostlearner, 
    resampling = rsmp("cv", folds = 5),
    measure = msr("classif.auc"), # ce
    search_space = search_space,
    terminator = trm("evals", n_evals = 10),    
    fitFunction = wfitXGBoost
  )
    
  result <- list(
    name = glue('xgboost max_depth {max_depth$param$lower}-{max_depth$param$upper}'),
    param = param,  
    fitFunction = wfitXGBoost
  )
  class(result) <- 'wmodel'
  return(result)
}


wfitXGBoost <- function(X, Y, param, w = NULL, verbose = F) {
  # todo change to formula
  if (verbose)
    if (!(missing(w) || is.null(w))) {
      cat('Fitting weighted XGBoost\n')
      print(summary(w*length(w)))
    }
  task <- mlr3preprocess(X, Y, w)
  
  instance = TuningInstanceSingleCrit$new(
    task = task,
    learner = param$learner,
    resampling = param$resampling,
    measure = param$measure,
    search_space = param$search_space,
    terminator = param$terminator
  )
  
  tuner = tnr("random_search")
  tuner$optimize(instance)
  param$learner$param_set$values = instance$result_learner_param_vals
  param$learner$predict_type = "prob"
  param$learner$train(task)
  result <- list(model=param$learner, predictFunction=wpredictXGBoost, importanceFunction=wimportantXGBoost)
  return(result)
}

wpredictXGBoost <- function(m, X) {
  # Notice X must be numeric
  p <- predict(m, as.data.frame(X), predict_type = "prob")[,2]
  return(p)
}


wimportantXGBoost <- function(m) {
  print(class(m$model$model))
  importantDf <- xgb.importance(model = m$model$model)
  important <- importantDf[['Feature']]
  return(important)
}


# TODO - this is copied and adapted from ml_utils.R remove duplication
mlr3preprocess <- function(X, Y, w) {
  # Notice X must be numeric
  dmat <- cbind(X, Y)
  d <- as.data.frame(dmat)
  d[['Y']] <- as.factor(d[['Y']])
  
  # TODO turning factors to numerics, understand this    
  indx <- sapply(d, is.factor)
  indx[['Y']] <- FALSE
  d[indx] <- lapply(d[indx], function(v) as.numeric(as.character(v)))
  
  if (missing(w) || is.null(w)) {
    task <- TaskClassif$new(id="XGB", backend = d,target = "Y")  # used to be makeClassifTask in mlr
  } else {
    d[['w']] <- w/sum(w)*length(w)
    task <- TaskClassif$new(id="XGB", backend = d,target = "Y")  # used to be makeClassifTask in mlr
    task$set_col_roles('w', roles = 'weight')
  }
  return (task)
}