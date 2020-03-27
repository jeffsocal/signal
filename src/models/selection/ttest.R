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

library(parallel)

parallelTTest <- function(v_features=c(),
                          c_predict='predict',
                          d_data=c(),
                          n_features=10,
                          returnTop=NULL,
                          returnType='table',
                          n_cores=1){

  if(is.null(returnTop))
    returnTop <- n_features

  returnTop <- min(returnTop, length(v_features))

  l_tt <- mclapply(1:length(v_features),
                   simpleTTest,
                   v_features=v_features,
                   c_predict=c_predict,
                   d_data=d_data,
                   mc.cores=n_cores)

  v_pvals <- unlist(l_tt)

  d_pvals <- data.frame(feature=v_features,pval=v_pvals)
  d_pvals <- d_pvals[order(d_pvals$pval),]
  d_pvals$p.adjust.bh <- p.adjust(d_pvals$pval, method='BH')

  if( grepl(returnTop,'significant') ){
    d_pvals <- d_pvals[d_pvals$p.adjust.bh <= 0.05,]
  } else if( is.numeric(returnTop) ) {
    d_pvals <- d_pvals[1:returnTop,]
  } else {
    stop(paste0("returnTop: '", returnTop, "' not valid, must be 'sig' or as.interger()"))
  }

  if( returnType == 'table')
    return(d_pvals)

  if( returnType == 'vector')
    return(as.character(d_pvals$feature))

}

simpleTTest <- function(x=1,
                        v_features,
                        c_feature=v_features[x],
                        c_predict,
                        d_data
) {

  predictors <- unique(d_data[,c_predict])
  v_class_a <- which(d_data[,c_predict] == predictors[1])
  v_class_b <- which(d_data[,c_predict] == predictors[2])

  v_var_a <- d_data[v_class_a, c_feature]
  v_var_b <- d_data[v_class_b, c_feature]

  result = tryCatch({
    l_ttest <- wilcox.test(v_var_a, v_var_b)
    this_pval <- l_ttest$p.val
    if( is.na(this_pval) | this_pval == 1) {
      this_pval <- 1
    }

    return(this_pval)
  }, warning = function(w) {
    return(1)
  }, error = function(e) {
    return(1)
  }, finally = {

  })


}
