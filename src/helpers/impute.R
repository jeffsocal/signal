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


static <- function(x, ...) {
  return(x)
}

impute <- function(m_data,
                   v_features,
                   method = median()){
  
  for ( fea in v_features ) {
    
    m_data[,fea] <- as.numeric(as.character(m_data[,fea]))
    
    this_zr <- which(m_data[,fea] == 0 | is.na(m_data[,fea]))
    
    if( length(this_zr) > 0 ) {
      this_fib <- length(this_zr)/dim(m_data)[1]
      # print.message('  impute zero values', paste0(fea, ": ",length(this_zr), " ", round(this_fib,3)*100, "%"))
      m_data[this_zr,fea] <- method(m_data[-this_zr,fea], na.rm=T)
    }
    
  }
  
  return(m_data)
  
}


