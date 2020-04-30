################################################################################
# Copyright (2016) SoCal Bioinformatics Inc. All rights reserved.
# This script is the confidential and proprietary product of SoCal
# Bioinformatics Inc. Any unauthorized reproduction or transfer of
# the contents herein is strictly prohibited.
#
################################################################################
# AUTH:     Jeff Jones | SoCal Bioinofrmatics Inc (SCBI)
# DATE:     2017.01.01
# OWNER:    SoCal Bioinofrmatics Inc
# PROJECT:  SCBI | classifiers
# DESC:     given a data matrix of disease/contol w/ measurements ~ classify
################################################################################

library(caret)
library(plyr)

signal_final_freq <- function(obj){
  
  d_data      <- obj$inputs$data
  c_predict   <- obj$inputs$predict
  v_features  <- obj$inputs$features
  
  f_select    <- obj$inputs$selection
  f_model     <- obj$inputs$model
  
  n_folds     <- obj$inputs$folds
  n_reps      <- obj$inputs$reps
  
  obj[['final']] <- list() 
  
  m_rf <- expand.grid(1:obj$inputs$folds,
                      1:obj$inputs$reps)
  i_rf <- 1:dim(m_rf)[1]
  
  l_freq <- list()
  for ( i in i_rf ){
    r <- m_rf[i,2]
    f <- m_rf[i,1]
    
    l_freq[[i]] <- obj$results[[r]][[f]]$selection
  }
  
  d_freq <- l_freq %>% 
    bind_rows() %>%
    group_by(feature) %>%
    summarise(
      num_rank = length(rank),
      frq_rank = length(rank) / length(l_freq),
      min_rank = min(rank)) %>%
    arrange(desc(frq_rank), min_rank)
  
  obj$final$features <- d_freq
  return(obj)
}

signal_final_model <- function(obj){
  
  obj <- signal_final_freq(obj)
  
  classifier_model  <- obj$inputs$model
  
  f_preprocess      <- obj$inputs$preprocess
  prd               <- obj$inputs$predictor
  dat               <- obj$inputs$data
  nfs               <- obj$inputs$nfeatures
  
  fea               <- obj$final$features$feature
  
  v_fold            <- obj$results[[1]]$foldAssigments[[1]]
  
  # d_train            <- dat[-v_fold,]
  # d_test             <- dat[v_fold,]
  
  d_train            <- d_test <- dat
  
  # preprocess on training data, apply to both
  f_preprocess      <- preProcess(d_train, c('center', 'scale'))
  d_tr              <- predict(f_preprocess, d_train)
  d_ts              <- predict(f_preprocess, d_test)
  
  obj$final$classifier <- list()
  for(nf in nfs){
    
    out_model <- classifier_model(v_features=fea[1:nf],
                                  c_predict=prd,
                                  d_train=d_tr,
                                  d_test=d_ts)
    
    cls <- list()
    cls$model <- out_model
    cls$cutoffs <- list()
    
    # OPTIMAL
    perf <- performance(out_model$pred, "cost")
    cls$cutoffs$optimal <- as.numeric(out_model$pred@cutoffs[[1]][which.min(perf@y.values[[1]])])
    
    # HIGHEST ACCURACY
    perf <- performance(out_model$pred, measure = "acc")
    ind <- which.max( slot(perf, "y.values")[[1]] )
    acc <- slot(perf, "y.values")[[1]][ind]
    cutoff <- slot(perf, "x.values")[[1]][ind]
    cls$cutoffs$accuracy <- as.numeric(cutoff)
    cls$cutoffs$max_accuracy <- as.numeric(acc)
    
    obj$final$classifier[[nf]] <- cls
    
  }
  
  return(obj)
  
}
