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
source("./src/plot_annotations.R")

plot_roc_individuals <- function(obj, annotate = T){
  
  sig_ann <- plot_annotations(obj)
  c_title <- paste(sig_ann$title, "Individual Feature Performace")
  c_annotation <- sig_ann$annotation
  c_title_sub <- sig_ann$title_sub
  
  d_data      <- obj$inputs$data
  c_predict   <- obj$inputs$predict
  v_features  <- obj$final$features
  
  m_roc <- as.matrix(d_data[,v_features])*-1
  m_obs <- matrix(as.numeric(d_data[,c_predict]), nrow=dim(m_roc)[1], ncol=dim(m_roc)[2])
  
  pred <- prediction(m_roc, m_obs)
  perf <- performance(pred, "tpr", "fpr")
  
  pauc <- performance(pred, "auc")
  
  obj_auc <- mean(unlist(pauc@y.values))
  obj_auc_sd <- sd(unlist(pauc@y.values))
  v_aucs <- signif(unlist(pauc@y.values),3)
  
  d_perf <- c()
  for( i in 1:length(v_aucs)){
    d_perf <- d_perf %>%
      bind_rows(
        tibble(
          fpr = perf@x.values[i] %>% unlist(),
          tpr = perf@y.values[i] %>% unlist(),
          feature = paste(v_features[i],
                          "AUC", v_aucs[i])
        )
      )
  }
  
  p_perf <- d_perf %>%
    ggplot(aes(fpr,tpr, color=feature)) +
    geom_abline(color='grey', linetype=2) +
    geom_step() +
    scale_x_continuous(breaks=c(0:5/10,.75,1), minor_breaks = NULL) +
    scale_y_continuous(breaks=c(5:10/10,.25,0), minor_breaks = NULL) +
    ggtitle(c_title, c_title_sub) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8)) +
    xlab(obj$performance$roc_perf@x.name) +
    ylab(obj$performance$roc_perf@y.name)
  
  if(annotate == T)
    p_perf <- p_perf +
    annotate('text', x=1,y=0, label=c_annotation,
             hjust=1, vjust=0, size=3)
  
  # plot(p_perf) 
  return(p_perf)
  
}
