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


# c_    character
# v_    vector
# d_    data.frame
# l_    list
# m_    matrix
# n_    numeric integer count
# i_    integer, index
# f_    float, function
# p_    probability
# === CODE =====================================================================

library(caret)
library(plyr)

# organize user imput
source("./src/signal_usr_input.R")

# AUTOLOAD HELPER FUNCTIONS
for( helper in list.files("./src/helpers/", full.names = T)){
  source(helper)
}

# FEATURE SELECTION                   bf enet ga gafeat mine rf ttest
source("./src/model_selection.R")

# CLASSIFIERS                         glm ksvm  rf  svm
source("./src/model_classification.R")


source("./src/signal_initialize.R")
source("./src/signal_build.R")
source("./src/signal_roc.R")
source("./src/signal_finalize.R")

# operations model
# ================
# REPS
# .   FOLDS
# .   .   train
# .   .   .   feature select
# .   .   .   .   MIC -> TTEST -> GA = FEATURES
# .   .   .   .   (whatever you can think of)
# .   .   .   model build
# .   .   .   .   SVM -> SUn = MODEL
# .   .   .   .   (whatever you can think of)
# .   .   test
# .   .   .   MODEL Predict
# .   .   .   ROC   ESTIMATE
# .   ROC-PERF      N-samples x 1
# ROC-PERF          N-samples x REPS


# data obj model
# =================
# obj                                  list
# ......inputs                            .list
# ..............dData                     ..matrix / data.frame
# ..............cPredictor                ..character (column in dData)
# ..............vfeatures                 ..vector
# ..............fSelection                ..list-function
# ..............fModel                    ..list-function
# ..............nFolds                    ..int (features)
# ..............nReps                     ..int (features)
# ..............nParallelCores            ..int (features)
# ......results                           .list
# ..............rep[i]                    ..list
# ...................foldAssignments      ...list
# ...................fold[i]              ...list
# .........................features       ....vector
# .........................trainAUC       ....numeric
# .........................testAUC        ....numeric
# .........................prediction     ....vector
# .........................observation    ....vector
# ......performance                       .ROCR.obj
# ......userInfo                          .list


signal <- function(
  mod_data = NULL,
  usr_project = 'SomeProject',
  usr_dataset = 'SomeData',
  usr_n_features = 3,
  usr_features = 'all',
  usr_cfv_split = NULL,
  usr_method = 'ttest',
  usr_model = 'svm05',
  usr_n_fold = 10,
  usr_n_reps = 10,
  usr_n_cores = 1,
  usr_verbose = F
){
  
  
  # override with command args if run from a script
  # however -- a params file should be constructed for each run
  cmd_args = commandArgs();
  for (arg in cmd_args){
    arg_value <- sub("--[a-z]*\\=", "", arg)
    if( grepl("--project", arg) ) usr_project <- arg_value
    if( grepl("--dataset", arg) ) usr_dataset <- arg_value
    if( grepl("--filter", arg) ) usr_features <- arg_value
    if( grepl("--method", arg) ) usr_method <- arg_value
    if( grepl("--model", arg) ) usr_model <- arg_value
    if( grepl("--cores", arg) ) usr_n_cores <- as.numeric(as.character(arg_value))
    if( grepl("--nfea", arg) ) usr_n_features <- as.character(arg_value)
    
    if( grepl("--nfold", arg) ) usr_n_fold <- as.character(arg_value)
    if( grepl("--nreps", arg) ) usr_n_reps <- as.character(arg_value)
    
    if( grepl("--split", arg) ) usr_cfv_split <- as.character(arg_value)
  } 
  
  user_info <- list()
  user_info['project'] <- usr_project
  user_info['dataset'] <- usr_dataset
  user_info['n_features'] <- usr_n_features
  user_info['features'] <- usr_features
  user_info['cfv_split'] <- usr_cfv_split
  user_info['method'] <- usr_method
  user_info['model'] <- usr_model
  user_info['n_fold'] <- usr_n_fold
  user_info['n_reps'] <- usr_n_reps
  user_info['n_cores'] <- usr_n_cores
  
  test_n <- c(2:20)
  if( grepl("^[0-9]+$", usr_n_features) ) {
    test_n <- as.numeric(usr_n_features)
  } else if( grepl(",", usr_n_features) ) {
    test_n <- as.numeric(unlist(strsplit(usr_n_features, ",")))
    test_n <- test_n[test_n <= length(fea)]
  } else if( grepl(":", usr_n_features) ) {
    test_n <- as.numeric(unlist(strsplit(usr_n_features, ":")))
    test_n <- test_n[1]:test_n[2]
  }
  
  if( usr_verbose == T){
    print_div()
    cat("USER SETTINGS\n")
    print_message('project', usr_project)
    print_message('dataset', usr_dataset)
    print_message('feature numbers', usr_n_features)
    print_message('     filter', usr_features)
    print_message('     selection', usr_method)
    print_message('     classification', usr_model)
    print_message('num. cfv folds/reps', paste0(usr_n_fold, "x", usr_n_reps))
    print_message('num. cpu cores', usr_n_cores)
    print_div()
    cat("DATA ATTRIBUTES\n")
    print_message('observations', dim(mod_data$data)[1])
    print_message('features', length(mod_data$features))
    print_div()
  }
  
  for ( nfeatures in test_n ){
    
    mat <- mod_data$data
    fea <- mod_data$features
    pre <- mod_data$preprocess
    prd <- mod_data$predict
    nfeatures <- min(nfeatures, length(fea))
    
    
    lm_start_time <- Sys.time()
    if( usr_verbose == T)
      cat("start", nfeatures, "features .. ")
    
    obj <- signal_initalize(d_data=mat,
                            c_predict=prd,
                            v_features=fea,
                            b_verbose=usr_verbose,  # <-- turn off for this type of iteration
                            n_reps=usr_n_reps,      # <-- adjust for testing 
                            n_folds=usr_n_fold,
                            n_features=nfeatures,
                            f_preprocess=pre,
                            n_cores=usr_n_cores,
                            l_usr_info=user_info,
                            c_split=usr_cfv_split
    )
    
    # set up the feature selection
    obj <- set_feature_selection(obj, usr_method)

    # set up the feature classification
    obj <- set_model_classification(obj, usr_model)
    
    # RUN the classification
    obj <- signal_build(obj)
    
    obj$timing$duration <- (Sys.time() - lm_start_time)
    
    auc <- round(median(obj$performance$roc_auc),3)
    
    if( usr_verbose == T)
      cat("\n  AUC:", auc, "..", obj$timing$duration, "sec \n")
    
    
    options(warn=-1)
    
    # SAVE the data
    saveResultObject(obj, 
                     usr_project, 
                     usr_dataset,
                     paste0(usr_features, 
                            "_cfv_", usr_n_fold, "x", usr_n_reps,
                            "_top_", nfeatures, 
                            "_", gsub(",", "-", usr_method), 
                            "_auc", auc))
    options(warn=0)
    
  }
  return(obj)
}


