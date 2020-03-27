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

# FUNCTIONS

library(e1071)
library(ROCR)
library(kernlab)

glmt <- function(v_features,
                 c_predict,
                 d_train,
                 d_test,
                 verbose=F,
                 type='glm',
                 ...){
  
  v_features <- v_features[!is.na(v_features)]
  
  out <- list(pred=NULL,
              auc=0,
              fpr=NULL,
              tpr=NULL,
              features=v_features,
              predictor=c_predict,
              prediction=NULL,
              observation=NULL)
  
  if(length(v_features) == 0)
    return(out)
  
  # if( grepl('svc', type ) ) {
  #   d_train[,c_predict] <- as.factor(d_train[,c_predict])
  #   d_test[,c_predict] <- as.factor(d_test[,c_predict])
  # } else {
  #   d_train[,c_predict] <- as.numeric(d_train[,c_predict] == 'disease')
  #   d_test[,c_predict] <- as.numeric(d_test[,c_predict] == 'disease')
  # }
  
  if( grepl('classification', type ) ) {
    d_train[,c_predict] <- as.factor(d_train[,c_predict])
    d_test[,c_predict] <- as.factor(d_test[,c_predict])
    # removed to insead use a predefined numeric class class
    # } else {
    #   d_train[,c_predict] <- as.numeric(d_train[,c_predict] == 'disease')
    #   d_test[,c_predict] <- as.numeric(d_test[,c_predict] == 'disease')
  }
  
  fmla <- paste(c_predict, "~ ", paste(v_features, collapse= "+"), sep="")
  fmla <- as.formula(fmla)
  
  # e1071::glm
  if ( type == 'glm' ){
    glm_model <- glm(fmla, data=d_train[,c(c_predict, v_features)])
    predict_model <- predict(glm_model, d_test[,c(c_predict, v_features)])
  }
  
  # kernlab::rvm
  if ( type == 'rvm' ){
    glm_model <- rvm(fmla, data=d_train[,c(c_predict, v_features)], verbosity=0)
    predict_model <- as.vector(predict(glm_model, d_test[,c(c_predict, v_features)]))
    names(predict_model) <- rownames(d_test)
  }
  
  if( is.null(predict_model) )
    return(out)
  
  #calculate the ROC AUC value
  pred     				<- prediction(predict_model, d_test[,c_predict])
  perf		 				<- performance(pred, "auc")
  
  out[['pred']] <- pred
  out[['perf']] <- perf
  out[['auc']] <- unlist(perf@y.values)
  
  perf    	 				<- performance(pred, "tpr", "fpr")
  out[['fpr']] <- unlist(perf@x.values)
  out[['tpr']] <- unlist(perf@y.values)
  
  out[['features']] <- v_features
  out[['predictor']] <- c_predict
  
  out[['classifier']] <- glm_model
  out[['prediction']] <- predict_model
  out[['observation']] <- d_test[,c_predict]
  
  if( verbose == T) {
    cat("--> GLM:", round(out[['auc']],3), ":")
    str(v_features)
  }
  
  return(out)
}
