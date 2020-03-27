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
library(ROCR)
library(randomForest)

rft <- function(v_features,
                 c_predict,
                 d_train,
                 d_test,
                 verbose=F,
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
  
  d_train[,c_predict] <- as.factor(d_train[,c_predict])
  d_test[,c_predict] <- as.factor(d_test[,c_predict])
  
  flma <- as.formula(paste0(c_predict, "~."))
  rf_model <- randomForest(flma, d_train[,c(c_predict, v_features)])
  
  predict_model <- predict(rf_model, d_test[,c(c_predict, v_features)])
  
  if( is.null(predict_model) )
    return(out)
  
  #calculate the ROC AUC value
  pred     				<- prediction(as.numeric(predict_model), d_test[,c_predict])
  perf		 				<- performance(pred, "auc")
  
  out[['pred']]        <- pred
  out[['perf']]        <- perf
  out[['auc']]         <- unlist(perf@y.values)
  
  perf    	 				   <- performance(pred, "tpr", "fpr")
  out[['fpr']]         <- unlist(perf@x.values)
  out[['tpr']]         <- unlist(perf@y.values)
  
  out[['features']]    <- v_features
  out[['predictor']]   <- c_predict
  
  out[['classifier']]  <- rf_model
  out[['prediction']]  <- predict_model
  out[['observation']] <- d_test[,c_predict]
  
  if( verbose == T) {
    cat("--> RF:", round(out[['auc']],3), ":")
    str(v_features)
  }
  
  return(out)
}

