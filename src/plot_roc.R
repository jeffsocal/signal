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

plot_roc <- function(obj, 
                     summary=c('reps','mean','smooth'),
                     annotate = T){
  
  sig_ann <- plot_annotations(obj)
  c_title <- paste(sig_ann$title,
                   obj$inputs$reps, "x", obj$inputs$folds, " CrossValidation")
  c_annotation <- sig_ann$annotation
  c_title_sub <- sig_ann$title_sub
  
  d_perf <- c()
  for( i in 1:length(obj$performance$roc_auc)){
    d_perf <- d_perf %>%
      bind_rows(
        tibble(
          fpr = obj$performance$roc_perf@x.values[i] %>% unlist(),
          tpr = obj$performance$roc_perf@y.values[i] %>% unlist(),
          rep = i
        )
      )
  }
  
  p_perf <- d_perf %>%
    arrange(fpr) %>%
    ggplot(aes(fpr,tpr)) +
    geom_abline(color='grey', linetype=2) +
    scale_x_continuous(breaks=c(0:5/10,.75,1), minor_breaks = NULL) +
    scale_y_continuous(breaks=c(5:10/10,.25,0), minor_breaks = NULL) +
    ggtitle(c_title, c_title_sub) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8)) +
    xlab(obj$performance$roc_perf@x.name[1]) +
    ylab(obj$performance$roc_perf@y.name[1])
  
  if('reps' %in% summary)
    p_perf <- p_perf + 
    geom_step(aes(group=rep),alpha=1/5)
  
  if('mean' %in% summary)
    p_perf <- p_perf + 
    geom_step(data = d_perf %>% group_by(fpr) %>% summarise(tpr=mean(tpr)),
              color='blue', size=1)
  
  if('smooth' %in% summary)
    p_perf <- p_perf + 
    stat_smooth(method='gam', color='blue', size=1, alpha=1/3)
  
  if(annotate == T)
    p_perf <- p_perf +
    annotate('text', x=1,y=0, label=c_annotation,
             hjust=1, vjust=0, size=3)
  
  # plot(p_perf) 
  return(p_perf)
  
}