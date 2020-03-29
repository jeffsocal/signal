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
  
  d_data[,c_predict] <- as.factor(d_data[,c_predict])
  
  d_freq <- c()
  # REPS
  for ( r in names(obj$results) ) {
    
    # FOLDS
    for ( f in names(obj$results[[r]][-1]) ) {
      
      d_freq <- rbind(d_freq,
                      data.frame(
                        feature=obj$results[[r]][[f]]$features,
                        replicate=r,
                        fold=f))
    }
  }
  
  d_f_rep <- ddply(d_freq, c('feature', 'replicate'), summarize, count=length(feature))
  d_f_all <- ddply(d_f_rep, c('feature'), summarize, count=median(count))
  d_f_rep$freq <- d_f_rep$count / n_reps
  
  d_f_all <- d_f_all[order(-d_f_all$count),]
  
  obj[['final']] <- list()
  obj$final$features <- as.character(d_f_all$feature[1:obj$inputs$nfeatures])
  
  return(obj)
}

signal_final_model <- function(obj){
  
  obj <- signal_final_freq(obj)
  
  classifier_model  <- obj$inputs$model
  
  f_preprocess      <- obj$inputs$preprocess
  prd               <- obj$inputs$predictor
  dat <- d_tr       <- obj$inputs$data
  fea               <- obj$final$features
  
  # d_tr              <- predict(f_preprocess, dat)
  
  out_model <- classifier_model(v_features=fea,
                                c_predict=prd,
                                d_train=d_tr,
                                d_test=d_tr)
  
  obj$final$classifier <- out_model$classifier
  obj$final[['cutoffs']] <- list()
  
  # OPTIMAL
  perf <- performance(out_model$pred, "cost")
  obj$final$cutoffs$optimal <- as.numeric(out_model$pred@cutoffs[[1]][which.min(perf@y.values[[1]])])
  
  # HIGHEST ACCURACY
  perf <- performance(out_model$pred, measure = "acc")
  ind <- which.max( slot(perf, "y.values")[[1]] )
  acc <- slot(perf, "y.values")[[1]][ind]
  cutoff <- slot(perf, "x.values")[[1]][ind]
  obj$final$cutoffs$accuracy <- as.numeric(cutoff)
  
  obj$final$max_accuracy <- as.numeric(acc)
  
  return(obj)
  
}
