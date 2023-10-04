library(magrittr)
library(tensorflow)
library(tfdatasets)
library(keras)
library(glue)


WeightedSumLogaritmic(keras$layers$Layer) %py_class% {
  initialize <- function(units = 1, input_dim = 32) {
    super$initialize()
    w_init <-  tf$constant_initializer(value = 0)
    self$w = self$add_weight(
      shape=shape(input_dim, units), initializer=w_init, trainable=T, constraint = constraint_nonneg()
    )
  }

  call <- function(inputs) {
    tf$matmul(inputs, tf$exp(self$w-tf$reduce_logsumexp(self$w)))
  }
}




WeightedSumId(keras$layers$Layer) %py_class% {
  initialize <- function(units = 1, input_dim = 32) {
    super$initialize()
    w_init <-  tf$constant_initializer(value = 1/input_dim)
    self$w = self$add_weight(
      shape=shape(input_dim, units), initializer=w_init, trainable=T, constraint = constraint_nonneg()
    )
  }

  call <- function(inputs) {
    tf$matmul(inputs, self$w)/tf$reduce_sum(self$w)
  }
}


optimizeKerasWeights <- function(wOptimizer, Z, mu) {
  # TODO - this is a duplication
  # Normalize if not all columns are binary
  isBinary <- apply(Z, 2,function(x) {all(x %in% c(0, 1))})
  if (!all(isBinary)) {
    normalized <- normalizeDataAndExpectations(Z, mu, wOptimizer$minSd)  # TODO
    rm(list=c('Z', 'mu'))
    gc()
    Z <- normalized$Z
    mu <- normalized$mu
  }
  Z <- as.matrix(Z)

  m <- ncol(Z)
  n <- nrow(Z)

  baseRates <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  if (wOptimizer$parameterization == 'id') {
    layer_weighted_sum <- create_layer_wrapper(WeightedSumId)
    rates <-  baseRates/1000
  } else {
    layer_weighted_sum <- create_layer_wrapper(WeightedSumLogaritmic)
    rates <- baseRates # *n*100
  }

  # callback <- keras$callbacks$EarlyStopping(monitor = 'loss', patience=10, baseline = -1, restore_best_weights = T)
  models <- vector(mode = 'list', length = length(rates))
  losses <- rep(NA, length(rates))
  for (i in 1:length(rates)) {
    if (wOptimizer$optimizer == 'adam')
      optimizer = optimizer_adam(learning_rate = rates[i])
    else
      optimizer = optimizer_sgd(learning_rate = rates[i], momentum = 0.9)

    model <- keras_model_sequential()
    model %>%
      layer_weighted_sum(units = 1, input_dim = n)

    model %>% compile(
      optimizer = optimizer,
      loss = 'mse'
    )

    model %>% fit(
      t(Z),
      mu,
      epochs = wOptimizer$nTuneIter,
      batch_size = m,
      verbose = 0
    )
    # print(model$history$history$loss) # [nIter]
    nIter <- length(model$history$history$loss)
    losses[i] <- model$history$history$loss[nIter]
    models[[i]] <- model
    cat(rates[i], nIter, losses[i], '\n')
  }
  j <- which.min(losses)
  cat('min' ,rates[j], losses[j], '\n')
  model <- models[[j]]
  # print(model$history$history$loss) # [nIter]

  model %>% fit(
    t(Z),
    mu,
    epochs = wOptimizer$nIter,
    batch_size = m,
    verbose = 0
    # callbacks = list(callback)
  )
  nIter <- length(model$history$history$loss)
  cat(nIter,'\n')
  print(model$history$history$loss[nIter]) #

  w_hat <- as.numeric(model$layers[[1]]$w)
  if (wOptimizer$parameterization != 'id')
    w_hat <- exp(w_hat)
  w_hat <- w_hat/sum(w_hat)

  # TODO duplicated
  muHat <- t(Z) %*% w_hat
  rr <- muHat[, 1] - mu
  err <- sum(rr**2)
  maxAbsErr <- max(abs(rr))

  if (err<wOptimizer$maxSuccessMSE)
    status = 'Success'
  else
    status = 'Not-converged'

  ParallelLogger::logInfo(glue('err {err} max abs err {maxAbsErr} status {status}'))
  return(list(model=model, w_hat = n * w_hat, err = err, maxAbsErr = maxAbsErr, status = status,
              totalIter=wOptimizer$nIter))
}


kerasWeightOptimizer <- function(
    optimizer, minSd=1e-4, nIter=5000, nTuneIter=50, outputDir=NULL, maxSuccessMSE=1e-5, parameterization = 'id') {
  l <- list(
    shortName = 'W-Keras',
    minSd = minSd,
    optimizer = optimizer,
    nIter = nIter,
    nTuneIter = nTuneIter,
    outputDir = outputDir,
    maxSuccessMSE = maxSuccessMSE,
    parameterization = parameterization,
    optimize = optimizeKerasWeights
  )
  class(l) <- 'kerasWeightOptimizer'
  return(l)
}
