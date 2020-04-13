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

# AUTOLOAD SELECTION FUNCTIONS
for( cf in list.files("./src/models/selection/", 
                      pattern = ".R", full.names = T)){
  source(cf)
}
ls_vars <- ls()

library(ROCR)
library(e1071)
library(stringr)
library(plyr)
library(doParallel)


feature_select.parallel <- function(i,
                                    obj,
                                    n_fun_cores = 1,
                                    ...){
  
  
  c_predict          <- obj$inputs$predict
  v_features         <- obj$inputs$features
  n_sel_fea          <- obj$inputs$nfeatures
  v_preprocess       <- obj$inputs$preprocess
  
  f_selects          <- obj$inputs$fselect
  f_selPars          <- obj$inputs$fselParallel
  n_fselTop          <- obj$inputs$fselTopN
  f_models           <- obj$inputs$model
  
  n_folds            <- obj$inputs$folds
  n_reps             <- obj$inputs$reps
  b_verbose          <- obj$inputs$verbose
  
  d_data             <- obj$inputs$data
  
  m_rf               <- expand.grid(1:obj$inputs$folds,1:obj$inputs$reps)
  rep                <- m_rf[i,2]
  fld                <- m_rf[i,1]
  
  v_fold             <- obj$results[[rep]]$foldAssigments[[fld]]
  
  d_train            <- d_data[-v_fold,]
  d_test             <- d_data[v_fold,]
  
  if(!'none' %in% v_preprocess){
    # preprocess on training data, apply to both
    f_preprocess       <- preProcess(d_train, v_preprocess)
    
    d_train            <- predict(f_preprocess, d_train)
    d_test             <- predict(f_preprocess, d_test)
  }
  
  # feature SELECTION METHODS: returns a tibble with columns; feature, rank
  v_fea_sel <- v_features
  for( i_fsel in 1:length(f_selects) ){
    
    f_select_name <- names(f_selects)[i_fsel]
    f_select <- f_selects[[f_select_name]]
    f_sel_par <- f_selPars[[f_select_name]]
    
    n_topFeatures <- max(n_sel_fea)*(length(f_selects)-(i_fsel-1))
    
    if( !is.na(n_fselTop[[f_select_name]]) )
      n_topFeatures <- n_fselTop[[f_select_name]]
    
    t_fea_sel <- f_select(v_features=v_fea_sel,
                          c_predict=c_predict,
                          d_data=d_train,
                          n_cores = f_sel_par)
    
    v_fea_sel <- t_fea_sel[1:n_topFeatures,]$feature %>% as.character()
    
  }
  return(t_fea_sel[1:max(n_sel_fea),])
}

feature_select <- function(obj,
                           fsel_parallel=feature_select.parallel,
                           ...) {
  
  m_rf <- expand.grid(1:obj$inputs$folds,1:obj$inputs$reps)
  i_rf <- 1:dim(m_rf)[1]
  
  # attempt at doParallel
  # cl <- makeCluster(2)
  # clusterExport(cl, ls_vars)
  # registerDoParallel(cl,library(thePackageYouUse))
  # 
  # l_features <- foreach(i=i_rf) %dopar% 
  #   fsel_parallel(i, obj)
  
  # one of the selection methods should advantage the cores
  n_sel_cores = obj$inputs$cores
  if(max(unlist(obj$inputs$fselParallel)) > 1)
    n_sel_cores = 1
  
  l_features <- mclapply(i_rf,
                         fsel_parallel,
                         obj=obj,
                         mc.cores=n_sel_cores)
  
  for ( i in i_rf ){
    rep <- m_rf[i,2]
    fld <- m_rf[i,1]
    obj$results[[rep]][[fld+1]]$selection <- l_features[[i]]
    obj$results[[rep]][[fld+1]]$classification <- list()
  }
  return(obj)
  
}

set_feature_selection <- function(obj,
                                  use_method="ttest:123"){
  
  # create a list of feature selection functions
  l_sel_func <- list(
    
    rf = function(...) rf_sel(
      ...),
    
    ttest = function(...) ttest_sel(
      ...),
    
    mine = function(...) mine_sel(
      ...),
    
    enet = function(...) enet_sel(
      ...),
    
    bf = function(...) bf_sel(
      returnType='vector',
      ...),
    
    ga = function(...) ga_sel(
      n_features=max(obj$inputs$nfeatures),
      n_max_generations=20,
      seed_method='rand',
      seed_n_comb=2,
      seed_n_limit=1e3,
      verbose=T,
      ...)
  )
  
  l_sel_returnTop <- list(
    rf = NA,
    ttest = NA,
    mine = NA,
    enet = NA,
    bf = NA,
    ga = NA
  )
  
  l_sel_parallel <- list(
    rf = 1,
    ttest = 1,
    mine = obj$inputs$cores,
    enet = 1,
    bf = obj$inputs$cores,
    ga = 1
  )
  
  for ( mth in unlist(strsplit(use_method, ",")) ){
    if( grepl(":", mth) ){
      mth <- unlist(strsplit(mth, ":"))
      l_sel_returnTop[[mth[1]]] <- as.numeric(as.character(mth[2]))
    }
  }
  fsel <- sub("\\:\\d+", "", unlist(strsplit(use_method, ",")))
  
  # feature SELECTION
  obj$inputs$fselect <- l_sel_func[fsel]
  obj$inputs$fselParallel <- l_sel_parallel[fsel]
  obj$inputs$fselTopN <- l_sel_returnTop[fsel]
  
  return(obj)
}
