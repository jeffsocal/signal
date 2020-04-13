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

signal_initalize <- function(d_data,
                             c_predict,
                             v_features,
                             n_folds=10,
                             n_reps=10,
                             n_cores=1,
                             f_preprocess=static, #impute
                             b_verbose=F,
                             n_features=5,
                             l_usr_info=list(),
                             c_split=NULL
){
  
  # if ( !is.object(f_preprocess) )
  #   stop("data preprocessing must be an object of type prediction")
  
  obj <- list(inputs=list(),
              results=list(),
              performance=list(),
              timing=list(start=Sys.time(),
                          end=Sys.time()),
              userInfo=l_usr_info)
  
  obj$inputs <- list(
    data=d_data,
    predictor=c_predict,
    features=v_features,
    fselect=list(),
    fselTopN=list(),
    fselParallel=list(),
    model=list(),
    folds=n_folds,
    reps=n_reps,
    cores=n_cores,
    preprocess=f_preprocess,
    verbose=b_verbose,
    nfeatures=n_features,
    split=c_split
  )
  
  # FOLD ASSIGNMENTS WITHIN REPS ===============================================
  l_folds <- list('foldAssigments'=c(0,0))
  for ( i in 1:n_folds ) {
    this_fld <- paste0("FOLD", str_pad(i, width=2, pad='0'))
    l_folds[[this_fld]] <- list()
  }
  
  
  for ( i in 1:n_reps ) {
    this_rep <- paste0("REP", str_pad(i, width=2, pad='0'))
    obj$results[[this_rep]] <- l_folds
    
    if(!is.null(c_split) & c_split %in% colnames(d_data) ){
      obj$inputs$split <- paste("xfolds split on id:", c_split)
      l_foldassign <- createFoldsFromID(d_data[,c_split],
                                        k = n_folds,
                                        list = T)
    } else {
      obj$inputs$split <- "xfolds split randomly"
      l_foldassign <- createFolds(d_data[,c_predict],
                                  k = n_folds,
                                  list = T)
    }
    obj$results[[this_rep]][['foldAssigments']] <- l_foldassign
  }
  # ============================================================================
  
  return(obj)
}
