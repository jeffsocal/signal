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

library(ROCR)
library(parallel)
library(e1071)
library(stringr)
library(plyr)

feature_select <- function(obj,
                           n_cores=detectCores()-1,
                           ...) {

  m_rf <- expand.grid(1:obj$inputs$folds,1:obj$inputs$reps)
  i_rf <- 1:dim(m_rf)[1]

  v_features <- mclapply(i_rf,
                         feature_select.parallel,
                         obj=obj,
                         mc.cores=n_cores)

  for ( i in i_rf ){
    rep <- m_rf[i,2]
    fld <- m_rf[i,1]
    pdm <- v_features[[i]]
    obj$results[[rep]][[fld+1]]$features <- as.character(unlist(pdm))

  }
  return(obj)

}

feature_select.parallel <- function(i,
                                    obj,
                                    ...){


  c_predict          <- obj$inputs$predict
  v_features         <- obj$inputs$features
  n_sel_fea          <- obj$inputs$nfeatures

  f_selects          <- obj$inputs$fselect

  n_fselTop          <- obj$inputs$fselTopN
  f_models           <- obj$inputs$model

  n_folds            <- obj$inputs$folds
  n_reps             <- obj$inputs$reps
  f_preprocess       <- obj$inputs$preprocess
  b_verbose          <- obj$inputs$verbose

  d_data             <- obj$inputs$data

  m_rf                  <- expand.grid(1:obj$inputs$folds,1:obj$inputs$reps)
  rep                   <- m_rf[i,2]
  fld                   <- m_rf[i,1]

  v_fold                <- obj$results[[rep]]$foldAssigments[[fld]]
  # d_train               <- predict(f_preprocess, newdata=d_data[-v_fold,])
  d_train               <- d_data[-v_fold,]


  # feature SELECTION METHODS
  v_fea_sel <- v_features
  for( i_fsel in 1:length(f_selects) ){

    f_select_name <- names(f_selects)[i_fsel]
    f_select <- f_selects[[f_select_name]]

    n_topFeatures <- n_sel_fea*(length(f_selects)-(i_fsel-1))

    if( !is.na(n_fselTop[[f_select_name]]) )
      n_topFeatures <- n_fselTop[[f_select_name]]

    v_fea_sel <- f_select(v_features=v_fea_sel,
                          c_predict=c_predict,
                          d_data=d_train,
                          n_features=n_topFeatures)

    # Remove weird manifestation of a rare bug
    v_fea_sel <- v_fea_sel[v_fea_sel != '(Other)']

    # if( b_verbose == T ){
    #   cat(paste0("   ",
    #              "REP", str_pad(rep, w=3),
    #              " FLD", str_pad(fld, w=3),
    #              " ", f_select_name, " ",
    #              str_pad(length(v_fea_sel), w=3),
    #              " >> "))
    #   str(v_fea_sel)
    # }

  }


  return(v_fea_sel)

}

# feature_select.validate <- function(obj,
#                                     d_val){
# 
# 
#   c_predict          <- obj$inputs$predict
#   v_features       <- obj$inputs$features
#   n_sel_fea          <- obj$inputs$nfeatures
# 
#   f_selects          <- obj$inputs$fselect
# 
#   n_fselTop          <- obj$inputs$fselTopN
#   f_models           <- obj$inputs$model
# 
#   n_folds            <- obj$inputs$folds
#   n_reps             <- obj$inputs$reps
#   preprocess         <- obj$inputs$preprocess
#   b_verbose          <- obj$inputs$verbose
# 
#   d_dis             <- obj$inputs$data
# 
#   d_train               <- d_dis
#   d_train[,v_features]  <- preprocess(d_train[,v_features])
# 
# 
#   # feature SELECTION METHODS
#   v_fea_sel <- v_features
#   for( i_fsel in 1:length(f_selects) ){
# 
#     f_select_name <- names(f_selects)[i_fsel]
#     f_select <- f_selects[[f_select_name]]
# 
#     n_topFeatures <- n_sel_fea*(length(f_selects)-(i_fsel-1))
# 
#     if( !is.na(n_fselTop[[f_select_name]]) )
#       n_topFeatures <- n_fselTop[[f_select_name]]
# 
#     v_fea_sel <- f_select(v_features=v_fea_sel,
#                           c_predict=c_predict,
#                           d_data=d_train,
#                           n_features=n_topFeatures)
# 
#     # Remove weird manifestation of a rare bug
#     v_fea_sel <- v_fea_sel[v_fea_sel != '(Other)']
# 
#   }
# 
#   return(v_fea_sel)
# 
# }



set_feature_selection <- function(obj,
                                  use_method="ttest:123"){

  # create a list of feature selection functions
  l_sel_func <- list(

    rf = function(...) rf(
      ...),

    rffs = function(...) rffs(
      ...),

    ttest = function(...) parallelTTest(
      returnType='vector',
      n_cores=1,
      ...),

    bf = function(...) parallelBF(
      returnType='vector',
      n_cores=1,
      ...),

    mine = function(...) parallelMINE(
      returnType='vector',
      n_cores=1,
      ...),

    enet = function(...) simpleENet(
      returnType='vector',
      n_cores=1,
      ...),

    ga = function(...) gafeat(
      n_max_generations=50,
      n_cores=1,
      seed_method='bf',
      seed_n_comb=2,
      seed_n_limit=1e3,
      verbose=T,
      ...)
  )

  l_sel_returnTop <- list(
    rf = NA,
    rffs = NA,
    ttest = NA,
    bf = NA,
    mine = NA,
    enet = NA,
    ga = NA
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
  obj$inputs$fselTopN <- l_sel_returnTop[fsel]

  return(obj)
}
