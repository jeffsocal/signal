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

library(parallel)

bf_sel <- function(v_features=c(),
                       c_predict='predict',
                       d_data=c(),
                       n_comb=2,
                       n_features=10,
                       returnTop=NULL,
                       returnType='table',
                       n_cores=1){

  if(is.null(returnTop))
    returnTop <- n_features

  n_population          <- choose(length(v_features), n_comb)
  returnTop             <- min(returnTop, length(v_features))
  m_combs               <- t(combn(v_features, n_comb))

  l_tt <- mclapply(1:n_population,
                   bf_svm,
                   m_combs=m_combs,
                   c_predict=c_predict,
                   d_data=d_data,
                   mc.cores=n_cores)
  v_auc <- unlist(l_tt)

  d_comb1 <- data.frame(feature=m_combs[,1],
                        auc=v_auc)
  d_comb2 <- data.frame(feature=m_combs[,2],
                        auc=v_auc)

  d_comb1 <- d_comb1[order(-d_comb1$auc),]
  d_comb2 <- d_comb2[order(-d_comb2$auc),]

  d_comb1$rank <- 1:dim(d_comb1)[1]
  d_comb2$rank <- 1:dim(d_comb2)[1]

  d_rank <- rbind(d_comb1, d_comb2)

  d_rank <- d_rank[order(d_rank$rank),]
  d_rank <- d_rank[-which(duplicated(d_rank$feature)),]
  rownames(d_rank) <- 1:dim(d_rank)[1]

  if( is.numeric(returnTop) ) {
    d_rank <- d_rank[1:returnTop,]
  } else {
    stop(paste0("returnTop: '", returnTop, "' not valid, must be 'sig' or as.interger()"))
  }

  if( returnType == 'table')
    return(d_rank)

  if( returnType == 'vector')
    return(as.character(d_rank$feature))

}

bf_svm <- function(x=1,
                      m_combs,
                      v_features=m_combs[x,],
                      c_predict,
                      d_data
) {

  source("./src/fclass/svm.R")

  out <- svmt(v_features=v_features,
              c_predict=c_predict,
              d_train=d_data,
              d_test=d_data,
              verbose=F)

  return(out[['auc']])

}
