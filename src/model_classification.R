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

# AUTOLOAD CLASSIFICATION FUNCTIONS
for( cf in list.files("./src/models/classification/", 
                      pattern = ".R", full.names = T)){
  source(cf)
}

library(ROCR)
library(parallel)
library(e1071)

model_classification <- function(obj,
                                 n_cores=detectCores()-1,
                                 ...) {
  
  # expand grid to reps|folds|features
  m_rf <- expand.grid(1:obj$inputs$folds,
                      1:obj$inputs$reps,
                      obj$inputs$nfeatures)
  
  i_rf <- 1:dim(m_rf)[1]
  
  out_models <- mclapply(i_rf,
                         model_classification.parallel,
                         obj=obj,
                         mc.cores=n_cores)
  
  for ( i in i_rf ){
    fld <- m_rf[i,1] 
    rep <- m_rf[i,2]
    nfs <- m_rf[i,3]
    pdm <- out_models[[i]]
    
    if( is.null(pdm$prediction) )
      return(obj)
    
    obj$results[[rep]][[fld+1]]$classification[[nfs]] <- pdm
    
  }
  
  rm(out_models)
  
  return(obj)
  
}

model_classification.parallel <- function(i,
                                          obj,
                                          ...){
  
  classifier_model  <- obj$inputs$model
  
  f_preprocess      <- obj$inputs$preprocess
  prd               <- obj$inputs$predictor
  dat               <- obj$inputs$data
  
  # expand grid to reps|folds|features
  m_rf <- expand.grid(1:obj$inputs$folds,
                      1:obj$inputs$reps,
                      obj$inputs$nfeatures)
  
  rep               <- m_rf[i,2]
  fld               <- m_rf[i,1]
  nfs               <- m_rf[i,3]
  
  v_fold            <- obj$results[[rep]]$foldAssigments[[fld]]
  fea               <- unlist(obj$results[[rep]][[fld+1]]$selection$feature)
  fea               <- fea[1:nfs]
  
  # preprocess on training data, apply to both
  f_preprocess      <- preProcess(dat[-v_fold,], c('center', 'scale'))
  
  d_tr              <- predict(f_preprocess, dat[-v_fold,])
  d_ts              <- predict(f_preprocess, dat[v_fold,])

  out_model <- classifier_model(v_features=fea,
                                c_predict=prd,
                                d_train=d_tr,
                                d_test=d_ts,
                                ...)
  
  return(out_model)
  
}


set_model_classification <- function(obj,
                                     use_model='ksvm61'){
  
  l_models <- c(svm_models(),
                ksvm_models(),
                list("knn" = function(...) knnc(...)),
                list("glm01" = function(...) glmt(...)),
                list("glm02" = function(...) glmt(type='rvm', ...)),
                list("rf01" = function(...) rft(type='rvm', ...))
  )
  
  obj$inputs$model <- l_models[[use_model]]
  
  return(obj)
}

list_models <- function(){
  l_models <- c(svm_models(),
                ksvm_models(),
                list("glm01" = function(...) glmt(...)),
                list("glm02" = function(...) glmt(type='rvm', ...)))
  
  return(as.character(names(l_models)))
}

