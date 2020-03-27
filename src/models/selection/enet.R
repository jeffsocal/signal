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

simpleENet <- function(v_features=c(),
                       c_predict='predict',
                       d_data=c(),
                       n_features=10,
                       returnTop=NULL,
                       returnType='table',
                       n_cores=1){

  if(is.null(returnTop))
    returnTop <- n_features

  returnTop             <- min(returnTop, length(v_features))

  d_data[,c_predict] <- as.numeric(d_data[,c_predict] == unique(d_data[,c_predict])[2] )
  obj_enet <<- enet(as.matrix(d_data[,v_features]), d_data[,c_predict])

  d_enet <- data.frame(feature=as.character(names(obj_enet$meanx)),
                       enet_score=obj_enet$normx)

  d_enet <- d_enet[order(-d_enet$enet_score),]


  if( is.numeric(returnTop) ) {
    d_enet <- d_enet[1:returnTop,]
  } else {
    stop(paste0("returnTop: '", returnTop, "' not valid, must be as.interger()"))
  }

  if( returnType == 'table')
    return(d_enet)

  if( returnType == 'vector')
    return(as.character(d_enet$feature))

}
