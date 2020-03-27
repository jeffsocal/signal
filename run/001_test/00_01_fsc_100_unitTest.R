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

rm(list=ls())
# import the code
source("./src/signal.R")

# define the data
r_data <- readRDS("./dat/000_unittest/univariate_100fea-100samples_2sig.R")

r_data$data[,"patient_integer"] <- as.numeric(r_data$data[,'patient_status'] == "disease")
r_data$data[,"patient_pid"] <- rep(1:(length(r_data$data[,"patient_integer"])/2),2)
fea <- r_data$features

wtf <- r_data$data[,c("patient_pid","patient_integer")]

# apply preprocessing 
modeling_data <- preprocess(r_data$data, fea, "patient_integer")

class_model <- signal(
  modeling_data,
  usr_project = 'unittest',
  usr_dataset = 'unittest_100x100',
  usr_n_features = '2',
  usr_features = 'all',
  usr_cfv_split = 'patient_pid',
  usr_method = 'ttest',
  usr_model = 'svm05',
  usr_n_fold = 10,
  usr_n_reps = 10,
  usr_n_cores = 1,
  usr_verbose = T
)
