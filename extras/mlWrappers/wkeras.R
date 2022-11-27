# See https://tensorflow.rstudio.com/

# reticulate with keras installation issue: see https://github.com/rstudio/keras/issues/615 answer by apantovic

library(reticulate)
library(tfdatasets)
library(keras)

wkeras <- function(mid_layer_params = NULL, suffix = NULL, epochs=NULL) {
  if (is.null(mid_layer_params))
    mid_layer_params <- list(activation = "relu", units = 32)
  if (is.null(epochs)) epochs <- 30
  param <- list(mid_layer_params = mid_layer_params, epochs = epochs)

  if (!missing(suffix) || !is.null(suffix))
    name <- glue('keras-{suffix}')
  else
    name <- 'keras'
  
  result <- list(
    name = name,
    param = param,
    fitFunction = wfitKeras
  )
  class(result) <- 'wmodel'
  return(result)
}


wfitKeras <- function(X, Y, param, w = NULL, verbose = F) {
  # todo change the API to formula?
  # prepare data
  if (verbose)
    if (!(missing(w) || is.null(w))) {
      cat('Fitting weighted keras\n')
      print(summary(w*length(w)))
    }
  d <- as_tibble(cbind(X, Y))
  xFeatures <- glue('X{1:ncol(X)}')
  names(d) <- c(xFeatures, 'Y')
  # spec
  spec <- feature_spec(d, Y ~ .) %>% 
    step_numeric_column(
      all_numeric(),
      normalizer_fn = scaler_standard()
    )
  spec_prep <- fit(spec)
  # str(spec_prep$dense_features())
  # Model
  input <- layer_input_from_dataset(d %>% select(-Y))
  output <- input %>% 
    layer_dense_features(dense_features(spec_prep)) %>% 
    (function(x) {do.call(layer_dense, c(list(x), param$mid_layer_params))}) %>%
    layer_dropout(0.2) %>%
    # layer_dense(units = 16, activation = 'relu') %>%
    # layer_dropout(0.1) %>%
    layer_dense(units = 1, activation = "sigmoid")  # TODO , kernel_regularizer=regularizer_l2(l=1e-1)
  model <- keras_model(input, output)
  # Training
  model %>% compile(
    loss = loss_binary_crossentropy, 
    optimizer = "adam", # optimizer_adam(learning_rate = 0.001),
    metrics = 'binary_accuracy',  # binary_crossentropy, binary_accuracy, metric_auc(from_logits = F)
    weighted_metrics =  metric_auc(from_logits = F)
  )

  fitVerbose <- 0
  if (!(missing(w) || is.null(w))) {
    history <- model %>% 
      fit(
        x = d %>% select(-Y),
        y = d$Y,
        sample_weight = w,
        epochs = param$epochs, 
        validation_split = 0.2,
        verbose = fitVerbose # , batch_size = 128
      )
  } else {
    history <- model %>% 
      fit(
        x = d %>% select(-Y),
        y = d$Y, 
        epochs = param$epochs, 
        validation_split = 0.2,
        verbose = fitVerbose # , batch_size = 128
      )
  }
  
  result <- list(model=model, predictFunction=wpredictKeras, importanceFunction=wimportantKeras, param=param)
  return(result)
}


wpredictKeras <- function(m, X) {
  d <- as_tibble(X)
  xFeatures <- glue('X{1:ncol(X)}')
  names(d) <- xFeatures
  p <- predict(m, d, verbose=0)
  p <- as.vector(p)
  return(p)
}


wimportantKeras <- function(m) {
  return(NULL)
}


wprintKeras <- function(m) {
  print(m$model)
}
