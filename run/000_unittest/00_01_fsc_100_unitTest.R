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
data <- readRDS("./dat/000_unittest/univariate_100fea-100samples_2sig.R")
data <- readRDS("./dat/000_unittest/binary-x_500fea-120samples.R")

dat <- data$data
tfea <- data$target_features
fea <- data$features
fea <- unique(c(tfea, sample(data$features, 10)))
dat <- dat[,c(fea, 'patient_status')]

# create a binary 0/1 patient status
dat <- cbind(dat, as.numeric(dat[,'patient_status'] == "disease"))
colnames(dat)[dim(dat)[2]] <- 'patient_integer'
dat <- cbind(dat, rep(1:(length(dat[,"patient_integer"])/2),2))
colnames(dat)[dim(dat)[2]] <- 'patient_pid'

# apply preprocessing 
modeling_data <- preprocess(dat, fea, "patient_integer")


class_model <- signal(
  modeling_data,
  usr_project = 'unittest',
  usr_dataset = 'unittest_100x100',
  usr_n_features = '2:3',
  usr_features = 'all',
  usr_cfv_split = 'patient_pid',
  usr_method = 'rf',
  usr_model = 'svm05',
  # usr_model = 'ksvm61',
  usr_n_fold = 5,
  usr_n_reps = 10,
  usr_n_cores = 1,
  usr_verbose = T
)


# usr_project = 'unittest'
# usr_dataset = 'unittest_100x100'
# usr_n_features = '2'
# usr_features = 'all'
# usr_cfv_split = 'patient_pid'
# usr_method = 'enet'
# usr_model = 'svm05'
# usr_n_fold = 10
# usr_n_reps = 10
# usr_n_cores = 1
# usr_verbose = T

