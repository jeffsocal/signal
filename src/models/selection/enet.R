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
# DESC:     function for incorporating elastic net regularization
#           NOTE: enet is only capable of finding models of linear combination
################################################################################


library(elasticnet)

enet_sel <- function(v_features=c(),
                       c_predict='predict',
                       d_data=c(),
                       returnType='table',
                       n_cores=1){


  d_data[,c_predict] <- as.numeric(d_data[,c_predict] == unique(d_data[,c_predict])[2] )
  obj_enet <- enet(as.matrix(d_data[,v_features]), d_data[,c_predict])

  out <- as_tibble(obj_enet[c('normx','meanx')]) %>% 
    mutate(feature = names(obj_enet$meanx)) %>% 
    arrange(desc(meanx)) %>% 
    mutate(rank = row_number())
  
  return(out)

}
