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

# c_    character
# v_    vector
# d_    data.frame
# l_    list
# m_    matrix
# n_    numeric integer count
# i_    integer, index
# f_    float, function
# p_    probability
# === CODE =====================================================================

library(caret)
library(plyr)

signal_build <- function(obj, 
                   f_models = c(list("svm01" = function(...) svmt(...)))
){

  # run feature selection
  obj <- feature_select(obj)

  # run classification
  obj <- model_classification(obj)

  #run ROC performance estimation
  obj <- signal_roc(obj)
  
  # build the final feature frequency
  obj <- signal_final_freq(obj)
  
  # build the final model
  obj <- signal_final_model(obj)
  
  obj$timing$end <- Sys.time()
  
  return(obj)
}


signal_functionControl <- function(obj, 
                                   purpose='model', 
                                   FUN, ...){
  
  obj$inputs[[purpose]] = function(x, ...) runif(x, ...)
  
  return(obj)
}

