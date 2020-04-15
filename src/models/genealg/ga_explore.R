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
# DESC:     
################################################################################

ga_best <- function(obj,
        gens=1:length(obj$generations),
        verbose=T){
  
  for ( g in gens ) {
    wtf <- which(obj$generations[[g]]$population[which(obj$generations[[g]]$fitness == max(obj$generations[[g]]$fitness)),] == 1)
    if( verbose == T ){
      cat("gen", str_pad(g, width=3), ": ")
      cat(round(max(obj$generations[[g]]$fitness),4))
      cat(" >> ")
      cat(paste(wtf, collapse=" "), "\n")
    }
  }
  
  out <- list(max=max(obj$generations[[g]]$fitness),
              min=min(obj$generations[[g]]$fitness),
              mean=mean(obj$generations[[g]]$fitness))
  
  return(out)
}


ga_features <- function(obj,
                        gen=length(obj$generations)
){
  
  this_gen <- obj$generations[[gen]]
  this_best <- which(this_gen$fitness == max(this_gen$fitness))
  this_feat <- which(this_gen$population[this_best,] == 1)
  v_features <- as.character(obj$inputs$features[this_feat])
  
  return(v_features)
}


