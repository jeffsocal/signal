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


preprocess <- function(m_data,
                       v_features,
                       c_predictor = 'patient_integer'){
  
  v_nonfeatures <- setdiff(colnames(m_data), v_features)
  
  for ( fea in v_features ) {
    
    # MAY produce WARNINGS "NAs introduced by coercion "
    capture.output(m_data[,fea] <- as.numeric(as.character(m_data[,fea])))
    
  }
  
  m_data <- m_data[,c(v_nonfeatures, v_features)]
  
  r_data <- list()
  r_data[['data']] <- m_data
  r_data[['preprocess']] <- preProcess(m_data, method=c("center", "scale"))
  r_data[['features']] <- v_features
  r_data[['predict']] <- c_predictor
  
  return(r_data)
}


