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

# GET: Command Arguments =======================================================

signal_usr_input <- function(
  usr_project = 'SomeProject',
  usr_dataset = 'SomeData',
  usr_n_features = 3,
  usr_features = 'all',
  usr_cfv_split = 'id',
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

  user_info['n_choose'] <- test_n  
  
  if( use_verbose == T){
    print.div()
    cat("USER SETTINGS\n")
    print_message('project', use_project)
    print_message('dataset', use_dataset)
    print_message('feature numbers', use_n_features)
    print_message('     filter', use_features)
    print_message('     selection', use_method)
    print_message('     classification', use_model)
    print_message('num. cfv folds/reps', paste0(use_n_fold, "x", use_n_reps))
    print_message('num. cpu cores', use_n_cores)
    print.div()
    cat("DATA ATTRIBUTES\n")
    print_message('observations', dim(r_data$data)[1])
    print_message('features', length(r_data$features))
    print_div()
  }
  
  return(user_info)
}