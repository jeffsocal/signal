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
# DESC:     install signal dependancies
################################################################################


install_signal <- function(){
  
  pkgs <- c('ggplot2','e1071','caret','plyr','tidy','ROCR','gridExtra',
            'parallel','stringr','kernlab','randomForest','elasticnet','GA',
            'doParallel','minerva','tidyverse','purrr','broom','pypr')
  
  install.packages(pkgs)

  }