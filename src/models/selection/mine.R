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
library(minerva)


mcl_mine <- function (i, dat, set_list) {
  a <- set_list[i,][1]
  b <- set_list[i,][2]
  x <- dat[,a]
  y <- dat[,b]
  w <- mine(x,y)

  linereg <- lm(x~y)
  w$SLOPE <- as.numeric(linereg$coef[2])
  w$FEATURE_A <- as.character(a)
  w$FEATURE_B <- as.character(b)

  return(w)
}

parallelMINE <- function(v_features=c(),
                         c_predict='predict',
                         d_data=c(),
                         n_features=10,
                         returnTop=NULL,
                         returnType='table',
                         n_cores=1,
                         c_metric="MIC"){

  if(is.null(returnTop))
    returnTop <- n_features

  returnTop             <- min(returnTop, length(v_features))

  set_list <- combn(v_features, 2)
  set_list <- t(set_list)
  i_set_list <- 1:dim(set_list)[1]

  d_mine <- mclapply(i_set_list,
                       mcl_mine,
                       dat=d_data,
                       set_list=set_list,
                       mc.cores=n_cores)

  d_mine <- do.call(rbind.data.frame, d_mine)
  d_mine <- d_mine[order(-d_mine[,c_metric]),]

  d_rank <- rbind(data.frame(feature=d_mine$FEATURE_A,
                             metric=d_mine[,c_metric],
                             rank=1:dim(d_mine)[1]),

                  data.frame(feature=d_mine$FEATURE_B,
                             metric=d_mine[,c_metric],
                             rank=1:dim(d_mine)[1])
                  )
  colnames(d_rank)[2] <- c_metric
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
