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

rf_sel <- function(v_features,
               c_predict,
               d_data,
               c_method="boot",
               f_function=rfFuncs,
               ...) {

  f_control <- rfeControl(functions = f_function,
                          method = c_method,
                          repeats = 3,
                          verbose = F,
                          allowParallel = FALSE)

  rfe_profile <- rfe(d_data[,v_features],
                     d_data[,c_predict],
                     sizes = length(v_features),
                     rfeControl = f_control)

  out <- as_tibble(rfe_profile$fit$importance) %>%
    mutate(feature = rownames(rfe_profile$fit$importance)) %>%
    select(
      feature,
      IncMSE = `%IncMSE`,
      IncNodePurity
      ) %>%
    mutate(rank = row_number())
  
  return(out)

}
