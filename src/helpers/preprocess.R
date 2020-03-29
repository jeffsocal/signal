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


# preprocess <- function(m_data,
#                        v_features,
#                        c_predictor = 'patient_integer'){
#   
#   r_data <- list()
#   r_data[['data']] <- m_data
#   r_data[['preprocess']] <- preProcess(m_data[,v_features], method=c("center", "scale"))
#   r_data[['features']] <- v_features
#   r_data[['predict']] <- c_predictor
#   
#   return(r_data)
# }

preprocess <- function(m_data,
                       v_features,
                       c_predictor = 'patient_integer'){
  
  m_data[,v_features] <- scale(m_data[,v_features])
  
  r_data <- list()
  r_data[['data']] <- m_data
  r_data[['preprocess']] <- NULL
  r_data[['features']] <- v_features
  r_data[['predict']] <- c_predictor
  
  return(r_data)
}

