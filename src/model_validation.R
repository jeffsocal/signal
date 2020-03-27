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
for( cf in list.files("./src/models/classification/", full.names = T)){
  source(cf)
}

library(ROCR)
library(parallel)
library(e1071)

model_validation <- function(obj,
                             fea,
                             d_val,
                             ...){
  
  classifier_model  <- obj$inputs$model
  
  preprocess        <- obj$inputs$preprocess
  prd               <- obj$inputs$predictor
  d_dis             <- obj$inputs$data
  
  d_tr              <- d_dis
  d_ts              <- d_val
  
  d_tr[,fea]        <- preprocess(d_tr[,fea])
  d_ts[,fea]        <- preprocess(d_ts[,fea])
  
  out_model <- classifier_model(v_features=fea,
                                c_predict=prd,
                                d_train=d_tr,
                                d_test=d_ts,
                                ...)
  
  return(out_model)
  
}
