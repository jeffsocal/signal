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
library(randomForest)

rf <- function(v_features,
               c_predict,
               d_data,
               n_features=10,
               c_method="boot",
               f_function=rfFuncs,
               ...) {

  # list of functions
  #     linear regression:        lmFuncs
  #     random forests:           rfFuncs
  #     naive Bayes:              nbFuncs
  #     bagged trees:             treebagFuncs
  #     caret's train function:   caretFuncs

  f_control <- rfeControl(functions = f_function,
                          method = c_method,
                          repeats = 10,
                          verbose = F,
                          allowParallel = FALSE)

  rfe_profile <- rfe(d_data[,v_features],
                     d_data[,c_predict],
                     sizes = n_features,
                     rfeControl = f_control)

  d_fea <- as.data.frame(
    summary(
      as.factor(
        rfe_profile$variables[rfe_profile$variables$Variables == n_features,]$var
        )
      )
    )
  colnames(d_fea)[1] <- 'count'
  v_fea_sel <- rownames(d_fea)[order(-d_fea$count)][1:n_features]

  return(as.character(v_fea_sel))

}

rffs <- function(v_features,
               c_predict,
               d_data,
               n_features=10,
               c_method="boot",
               f_function=rfFuncs,
               ...) {

  # list of functions
  #     linear regression:        lmFuncs
  #     random forests:           rfFuncs
  #     naive Bayes:              nbFuncs
  #     bagged trees:             treebagFuncs
  #     caret's train function:   caretFuncs

  f_control <- rfeControl(functions = f_function,
                          method = c_method,
                          repeats = 10,
                          verbose = F,
                          allowParallel = FALSE)

  rfe_profile <- rfe(d_data[,v_features],
                     d_data[,c_predict],
                     sizes = 1:20,
                     rfeControl = f_control)


  return(as.character(rfe_profile$optVariables))

}
