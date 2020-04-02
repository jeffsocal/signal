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

library(e1071)
library(ROCR)

svm_cls <- function(v_features,
                 c_predict,
                 d_train,
                 d_test,
                 verbose=F,
                 type='eps-regression',
                 kernel='radial'){
  
  out <- list()
  v_features <- v_features[!is.na(v_features)]
  if(length(v_features) == 0){
    out[['auc']] <- 0
    return(out)
  }
  
  if( grepl('classification', type ) ) {
    d_train[,c_predict] <- as.factor(d_train[,c_predict])
    d_test[,c_predict] <- as.factor(d_test[,c_predict])
  }
  
  svm_model <- svm(d_train[,v_features], 
                   d_train[,c_predict], 
                   type=type, 
                   kernel=kernel, 
                   scale=FALSE)
  
  predict_model <- predict(svm_model, 
                           newdata=d_test[,v_features])
  
  #calculate the ROC AUC value
  pred <- prediction(predict_model, d_test[,c_predict])
  perf <- performance(pred, "auc")
  
  out[['pred']] <- pred
  out[['perf']] <- perf
  out[['auc']] <- unlist(perf@y.values)
  
  perf    	 				<- performance(pred, "tpr", "fpr")
  out[['fpr']] <- unlist(perf@x.values)
  out[['tpr']] <- unlist(perf@y.values)
  
  out[['features']] <- v_features
  out[['predictor']] <- c_predict
  
  out[['classifier']] <- svm_model
  out[['prediction']] <- predict_model
  out[['observation']] <- d_test[,c_predict]
  
  # if( verbose == T) {
  #   cat("--> SVM:", round(out[['auc']],3), ":")
  #   str(v_features)
  # }
  
  return(out)
}

# svmtcv <- function(..., f=svmt,
#                    d_data,
#                    c_predict,
#                    n_fold=10,
#                    n_times=1) {
#   
#   f_train <- createFolds(d_data[,c_predict], k = n_fold, list = T)
#   
#   v_auc <- list()
#   for ( i in 1:length(f_train) ){
#     
#     this_train <- d_data[-f_train[[i]],]
#     this_test <- d_data[f_train[[i]],]
#     
#     this_model <- f(...,
#                     c_predict=c_predict,
#                     d_train=this_train,
#                     d_test=this_test)
#     
#     v_auc[[i]] <- this_model
#     
#   }
#   return(v_auc)
# }
# 
# 
# 
# svmtcv.stats <- function(svm_cv, f=mean){
#   this_auc <- c()
#   for( i in 1:length(svm_cv)){
#     this_auc <- append(this_auc, svm_cv[[i]]$auc)
#   }
#   return(f(this_auc))
# }

svm_models <- function(){
  
  type <- c('C-classification', 'eps-regression', 'nu-regression')
  kernel <- c('linear','radial','polynomial')
  
  l_models <- list()
  l_models[['svm01']] <- function(...) svm_cls(..., type='eps-regression', kernel='linear')
  l_models[['svm02']] <- function(...) svm_cls(..., type='eps-regression', kernel='radial')
  l_models[['svm03']] <- function(...) svm_cls(..., type='eps-regression', kernel='polynomial')
  l_models[['svm04']] <- function(...) svm_cls(..., type='nu-regression', kernel='linear')
  l_models[['svm05']] <- function(...) svm_cls(..., type='nu-regression', kernel='radial')
  l_models[['svm06']] <- function(...) svm_cls(..., type='nu-regression', kernel='polynomial')
  #   l_models[['svm07']] <- function(...) svmt(..., type='C-classification', kernel='linear')
  #   l_models[['svm08']] <- function(...) svmt(..., type='C-classification', kernel='radial')
  #   l_models[['svm09']] <- function(...) svmt(..., type='C-classification', kernel='polynomial')
  
  
  return(l_models)
}
