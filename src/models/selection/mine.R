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
  
  w <- w %>% as_tibble() %>%
    arrange(desc(MIC), desc(MAS)) %>%
    head(100)
  
  return(w)
}

mine_sel <- function(v_features=c(),
                     c_predict='predict',
                     d_data=c(),
                     n_cores=1){
  
  set_list <- combn(v_features, 2)
  set_list <- t(set_list)
  i_set_list <- 1:dim(set_list)[1]
  
  d_mine <- mclapply(i_set_list,
                     mcl_mine,
                     dat=d_data,
                     set_list=set_list,
                     mc.cores = n_cores
  )
  
  d_mine <- d_mine %>% 
    bind_rows() %>% 
    arrange(desc(MIC)) %>%
    mutate(rank = row_number())  %>%
    gather(FEATURE_A, FEATURE_B, key='AB', value='feature') %>%
    group_by(feature) %>%
    summarise(pair_rank = min(rank),
              mic = max(MIC),
              mas = max(MAS),
              mev = max(MEV),
              mcn = max(MCN),
              mic_r2 = max(`MIC-R2`)) %>%
    arrange(desc(mic), desc(mas)) %>%
    mutate(rank = row_number())
  
  return(d_mine)
}
