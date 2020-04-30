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
library(ROCR)
library(tidyverse)

signal_roc <- function(obj){
  
  d_data              <- obj$inputs$data
  c_predict           <- obj$inputs$predict
  v_components        <- obj$inputs$components
  
  f_select            <- obj$inputs$selection
  f_model             <- obj$inputs$model
  
  n_folds             <- obj$inputs$folds
  n_reps              <- obj$inputs$reps
  n_feas              <- obj$inputs$nfeatures
  
  d_data[,c_predict]  <- as.factor(d_data[,c_predict])
  
  for(n_fea in n_feas){
    d_pred <- list()
    # REPS
    for ( r in names(obj$results) ) {
      v_smpl <- c()
      v_pred <- c()
      
      # FOLDS
      for ( f in names(obj$results[[r]][-1]) ) {
        this_pred <- obj$results[[r]][[f]]$classification[[n_fea]]$prediction
        
        smpl <- which(rownames(d_data) %in% (names(unlist(this_pred))))
        pred <- as.vector(unlist(this_pred))
        v_smpl <- c(v_smpl, as.numeric(as.character(smpl)) )
        v_pred <- c(v_pred, pred )
        
      }
      
      d_pred[[r]] <- data.frame(
        sample = v_smpl, 
        pred = v_pred,
        rep = r) %>% 
        mutate(rep = as.character(rep))
      
    }
    
    d_pred_names <- names(d_pred)
    
    d_pred <- d_pred %>% 
      bind_rows() %>%
      mutate(row_id = dplyr::row_number())
    
    d_pred <- d_pred %>%
      mutate(obs = d_data[d_pred$sample,c_predict]) %>%
      spread(key='rep', value='pred') %>%
      as.data.frame()
    
    m_roc <- as.matrix(d_pred[,d_pred_names])
    suppressWarnings(storage.mode(m_roc) <- "numeric")
    m_obs <- matrix(as.numeric(d_pred[,'obs']), nrow=dim(m_roc)[1], ncol=dim(m_roc)[2])
    
    pred <- prediction(m_roc, m_obs)
    perf <- performance(pred, "tpr", "fpr")
    
    pauc <- performance(pred, "auc")
    
    obj_auc <- mean(unlist(pauc@y.values))
    obj_auc_sd <- sd(unlist(pauc@y.values))
    
    obj[['performance']][[n_fea]] <- list()
    obj[['performance']][[n_fea]]$matrix_predictions <- d_pred
    obj[['performance']][[n_fea]]$roc_pred <- pred
    obj[['performance']][[n_fea]]$roc_perf <- perf
    obj[['performance']][[n_fea]]$roc_auc <- unlist(pauc@y.values)
  }
  return(obj)
}

