library(tidymodels)
library(tidyverse)



setglmnet <- function(mixture=0.5) {
  return(list(mixture=mixture))
}


fitglmnet <- function(trainer, d, w) {
  # model <- fit(trainer, Y ~ ., data = d) 

  trainer <- 
    logistic_reg(penalty = tune(), mixture = 0.5) %>% 
    set_engine("glmnet", alpha=0.1, weights = w)
  
  folds <- vfold_cv(d, v=5)
  
  my_grid <- tibble(penalty = 10^seq(-4, -2, length.out = 20))
  rec <- recipe(Y ~ ., data = d)
  
  wf <- workflow() %>%
    add_model(trainer) %>%
    add_recipe(rec)
  
  my_res <- wf %>% 
    tune_grid(resamples = folds,
              grid = my_grid,
              control = control_grid(verbose = FALSE, save_pred = TRUE),
              metrics = metric_set(accuracy))
  
  best_mod <- my_res %>% select_best("accuracy")
  final_fitted <- 
    finalize_workflow(wf, best_mod) %>%
    fit(data = d)
  return(final_fitted)
}


wglmnetTest <- function() {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  source(file.path(script_dir,'../simulations/riskModel2.R'))
  hyperParams <- getDefaultHyperParams()
  params <- modelParams(hyperParams)
  d <- sampleModel(params, 1000)
  d[['Y']] <- as.factor(d[['Y']])
  
  trainer <- setglmnet()
  print(trainer$mixture)
  w <- 1:nrow(d)
  w <- (w/sum(w))*nrow(d)
  m <- fitglmnet(trainer, d, w)
  
  print(tidy(m))
  pt <- predict(m, d[,1:3], type='prob')
  p <- pull(pt, var='.pred_1')
  print(head(p))
}