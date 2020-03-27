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

plot_freqrank <- function(obj, annotate = F){
  
  
  d_data      <- obj$inputs$data
  c_predict   <- obj$inputs$predict
  v_features  <- obj$inputs$features
  
  f_select    <- obj$inputs$selection
  f_model     <- obj$inputs$model
  
  n_folds     <- obj$inputs$folds
  n_reps      <- obj$inputs$reps
  
  d_data[,c_predict] <- as.factor(d_data[,c_predict])
  
  d_freq <- c()
  # REPS
  for ( r in names(obj$results) ) {
    
    # FOLDS
    for ( f in names(obj$results[[r]][-1]) ) {
      
      d_freq <- rbind(d_freq,
                      data.frame(
                        feature=obj$results[[r]][[f]]$features,
                        replicate=r,
                        fold=f))
    }
  }
  
  d_f_rep <- ddply(d_freq, c('feature', 'replicate'), summarize, count=length(feature))
  d_f_all <- ddply(d_f_rep, c('feature'), summarize, count=median(count))
  d_f_rep$freq <- d_f_rep$count / n_reps
  
  d_f_all <- d_f_all[order(-d_f_all$count),]
  d_f_all$feature <- factor(d_f_all$feature, levels=d_f_all$feature)
  d_f_rep$feature <- factor(d_f_rep$feature, levels=d_f_all$feature)

  sig_ann <- plot_annotations(obj)
  c_title <- paste(sig_ann$title, "Frequency Rank Plot")
  c_annotation <- sig_ann$annotation
  c_title_sub <- paste("Distribution Among", n_reps, "Replicates")
  
  p_rep <- ggplot(d_f_rep, aes(feature, freq)) +
    geom_boxplot(alpha=1/3) +
    stat_summary(fun=median, geom="line", aes(group=1), color='blue') +
    stat_summary(fun=median, geom="point", color='blue') +
    ggtitle(c_title, c_title_sub) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8)) +
    axis_x_rotate(60) + xlab("feature") + ylab("Frequency")
  
  if(annotate == T)
    p_rep <- p_rep +
    annotate('text', x=1,y=0, label=c_annotation,
             hjust=1, vjust=0, size=3)
  
  out <- list(rep=d_f_rep,
              all=d_f_all,
              p_rep=p_rep)
  
  return(out)
  
}
