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


data_prep <- function(m_data,
                      v_features,
                      c_predictor = 'patient_integer',
                      v_preprocess = c('center', 'scale')){
  
  r_data <- list()
  r_data[['data']] <- m_data
  r_data[['preprocess']] <- v_preprocess
  r_data[['features']] <- v_features
  r_data[['predict']] <- c_predictor
  
  return(r_data)
}

