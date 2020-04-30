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
library(pROC)
source("./src/plot_annotations.R")

plot_roc_final <- function(obj, 
                     features = 2,
                     summary=c('reps','mean'),
                     conf_level = .95,
                     annotate = F){
  
  sig_ann <- plot_annotations(obj, features)
  c_title <- paste0("Bootstrap of Final Model | ", sig_ann$title)
  c_annotation <- sig_ann$annotation
  c_title_sub <- ''#sig_ann$title_sub
  
  # estimate the confidence intervals  
  proc_obj <- roc(obj$final$classifier[[features]]$model$observation, 
                  obj$final$classifier[[features]]$model$prediction, 
                  ci=TRUE, plot=FALSE)
  
  auc_ci <- proc_obj$ci
  
  proc_obj_ci <- ci.se(proc_obj, specificities=seq(0, 1, l=25))
  df_obj_ci <- data.frame(fpr = as.numeric(rownames(proc_obj_ci)),
                          tpr = 0,
                          tpr_lower = proc_obj_ci[, 1],
                          tpr_upper = proc_obj_ci[, 3])
  
  c_title_sub <- paste0(c_title_sub, 
                        "AUC mean: ", signif(auc_ci[2],3), "  ",
                        signif(conf_level * 100,3),"% CI: ",
                        signif(auc_ci[1],3)," - ", signif(auc_ci[3],3)
  )
  
  p_perf <- proc_obj %>%
    ggroc() +
    geom_abline(slope=1, intercept = 1, color='grey', linetype=2) +
    ggtitle(c_title, c_title_sub) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8)) +
    xlab(paste("1 -", obj$performance[[features]]$roc_perf@x.name[1])) +
    ylab(obj$performance[[features]]$roc_perf@y.name[1]) +
    geom_ribbon(data = df_obj_ci, 
                aes(x = fpr,
                    ymin = tpr_lower, 
                    ymax = tpr_upper), fill = "steelblue", alpha= 0.2)
  
  if(annotate == T)
    p_perf <- p_perf +
    annotate('text', x=0, y=0, label=c_annotation,
             hjust=1, vjust=0, size=3)
  
  # plot(p_perf) 
  return(p_perf)
  
}