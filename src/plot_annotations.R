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

plot_annotations <- function(obj){
  out <- list()
  out$title <- paste(obj$userInfo$project)
  
  out$annotation <- ''
  for ( i in names(obj$userInfo) ){
    out$annotation <- paste0(out$annotation, " ", i, ": ",
                           paste(unlist(obj$userInfo[[i]]), collapse=" "),
                           "\n")
  }
  out$annotation <- sub("\n$", "", out$annotation)
  
  out$title_sub <- paste(
    "AUC ",
    "MEDIAN:", round(median(obj$performance$roc_auc),3),
    " MEAN:", round(mean(obj$performance$roc_auc),3),
    "-/+ ", round(sd(obj$performance$roc_auc),3)
  )
  
  return(out)
}
