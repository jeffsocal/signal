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
source("./src/plot_annotations.R")

plot_roc <- function(obj, 
                     features = 2,
                     summary=c('reps','mean'),
                     conf_level = 0.95,
                     annotate = F){
  
  conf_level_upper <- 1 - (1-conf_level) / 2
  conf_level_lower <- (1-conf_level) / 2
  
  sig_ann <- plot_annotations(obj, features, conf_level)
  c_title <- paste0("CrossValidation ", obj$inputs$reps, "x", obj$inputs$folds,
                    " | ", sig_ann$title)
  c_annotation <- sig_ann$annotation
  c_title_sub <- sig_ann$title_sub
  
  d_perf <- c()
  for( i in 1:length(obj$performance[[features]]$roc_auc)){
    d_perf <- d_perf %>%
      bind_rows(
        tibble(
          fpr = obj$performance[[features]]$roc_perf@x.values[i] %>% unlist(),
          tpr = obj$performance[[features]]$roc_perf@y.values[i] %>% unlist(),
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
    xlab(obj$performance[[features]]$roc_perf@x.name[1]) +
    ylab(obj$performance[[features]]$roc_perf@y.name[1])
  
  if('reps' %in% summary)
    p_perf <- p_perf + 
    geom_step(aes(group=rep), color='grey', 
              alpha=1/5)
  
  if('mean' %in% summary){
    
    tmp_perf <- d_perf %>% 
      # mutate(fpr = round(log10(fpr),1)) %>%
      # mutate(fpr = 10^fpr) %>%
      group_by(fpr) %>% 
      summarise(tpr_min=quantile(tpr,conf_level_lower)[1],
                tpr_max=quantile(tpr,conf_level_upper)[1],
                tpr = median(tpr)) %>%
      ungroup()
    
    p_perf <- p_perf + 
      geom_ribbon(data = tmp_perf,
                  aes(x=fpr, ymin=tpr_min, ymax=tpr_max), 
                  fill='steelblue', alpha=1/5) + 
      geom_step(data = tmp_perf,
                  color='red', alpha=1/2)
  }
  
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