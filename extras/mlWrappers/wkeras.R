# See https://tensorflow.rstudio.com/

# reticulate with keras installation issue: see https://github.com/rstudio/keras/issues/615 answer by apantovic

library(R6)
library(reticulate)
library(tfdatasets)
library(keras)

wkeras <- R6Class(
  "wkeras",

  public = list(

    epochs = NA,
    step_epochs = NA,
    name = 'keras',
    w = 0,
    p = NA,
    model = NULL,
    batch_size = NA,
    lambda = NA,
    nhidden = 0,

    #' @param epoches number of epochs in the initial training phase
    #' @param step_epoches number of epoths in the text training step
    initialize = function(epochs = 30, step_epochs=5, batch_size=128, lambda=0.01, nhidden=0) {
      stopifnot(is.numeric(epochs), length(epochs) == 1)
      self$epochs <- epochs
      stopifnot(is.numeric(step_epochs), length(step_epochs) == 1)
      self$step_epochs <- step_epochs
      self$batch_size <- batch_size
      self$lambda <- lambda
      self$nhidden <- nhidden
      self$name <- glue('wkeras-{nhidden}')
    },

    #' @param param for backward compatability
    fitFunction = function(X, Y, param, w = NULL, verbose = F) {

      # prepare data
      if (verbose)
        if (!(missing(w) || is.null(w))) {
          cat('Fitting weighted keras\n')
          print(summary(w*length(w)))
        }

      epochs <- self$step_epochs
      # init model
      if (is.na(self$p)) {
        #
        self$p <- ncol(X)
        ParallelLogger::logInfo(glue('Intializing keras p = {self$p}, lambda = {self$lambda}'))
        # init model
        if (self$nhidden == 0)
          self$model <- keras_model_sequential() %>%
            layer_dense(units = 1,
                        activation = "sigmoid",
                        kernel_regularizer = regularizer_l1_l2(l1 = self$lambda/2, l2 = self$lambda/2),
                        input_shape = shape(self$p))
        else
          self$model <- keras_model_sequential() %>%
          layer_dense(units = 64,
                      activation = "relu",
                      kernel_regularizer = regularizer_l1_l2(l1 = self$lambda/2, l2 = self$lambda/2),
                      input_shape = shape(self$p)) %>%
          layer_dense(units = 1,
                      activation = "sigmoid",
                      kernel_regularizer = regularizer_l1_l2(l1 = self$lambda/2, l2 = self$lambda/2),
                      input_shape = shape(64))

        # Training
        self$model %>% compile(
          loss = "binary_crossentropy", # loss_binary_crossentropy,
          optimizer = "adam", # optimizer_adam(learning_rate = 0.001),
          metrics = c("accuracy"), # , 'binary_accuracy',  binary_crossentropy, binary_accuracy, metric_auc(from_logits = F)
          weighted_metrics =  c("accuracy")
        )
        epochs <- self$epochs
      }

      callbacks <- list(
        callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 5, verbose = 1, mode = "auto"))
      # callback_reduce_lr_on_plateau()

      fitVerbose <- 0
      if (!(missing(w) || is.null(w))) {
        history <- self$model %>%
          fit(
            x = as_tensor(X),
            y = as_tensor(Y),
            sample_weight = w/sum(w)*length(w),
            epochs = epochs,
            validation_split = 0.2,
            verbose = fitVerbose,
            batch_size = self$batch_size,
            callbacks = callbacks
          )
      } else {
        history <- self$model %>%
          fit(
            x = as_tensor(X),
            y = as_tensor(Y),
            epochs = epochs,
            validation_split = 0.2,
            verbose = fitVerbose,
            batch_size = self$batch_size,
            callbacks = callbacks
          )
      }
      return(self)
    },

    predictFunction = function(m=NULL, X) {  # TODO - get rid of the first argument, turn all classes to R6
      d <- as_tibble(X)
      xFeatures <- glue('X{1:ncol(X)}')
      names(d) <- xFeatures
      p <- predict(self$model, as_tensor(d), verbose=0)
      p <- as.vector(p)
      return(p)
    },


    importanceFunction = function(m=NULL) {
      return(NULL)
    },

    printFucntion = function(m=NULL) {
      print(self$model)
    }

  )
)
