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
# DESC:     plot various 
################################################################################

library(ggplot2)
library(tidyverse)
library(cvAUC)

plot_annotations <- function(obj, features = 2, conf_level = 0.9){
  out <- list()
  out$title <- paste(obj$userInfo$project)
  
  out$annotation <- ''
  for ( i in names(obj$userInfo) ){
    out$annotation <- paste0(out$annotation, " ", i, ": ",
                           paste(unlist(obj$userInfo[[i]]), collapse=" "),
                           "\n")
  }
  out$annotation <- sub("\n$", "", out$annotation)
  
  
  auc_ci <- ci.cvAUC(obj$performance[[features]]$roc_pred@predictions, 
                     obj$performance[[features]]$roc_pred@labels,
                     confidence = conf_level)
  
  out$title_sub <- paste0(
    "AUC ", signif(auc_ci$cvAUC,3), "  ",
    signif(conf_level* 100, 2), "% CI:", 
    signif(auc_ci$ci[1],3), " - ", signif(auc_ci$ci[2],3)
  )
  
  return(out)
}
