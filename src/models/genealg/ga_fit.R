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

# public (parallelize this function)
ga_fit <- function(obj, 
                   ...){
  
  # use the max_features to modify the ranking fitness of each feature
  n_max_features      <- obj$settings$max_features
  n_fit_penalty       <- obj$settings$fit_penalty
  c_fit_metric        <- obj$settings$fit_metric
  n_generations       <- obj$settings$generations
  n_cores             <- obj$settings$cores
  v_features          <- obj$inputs$features
  c_predict           <- obj$inputs$predict
  d_data              <- obj$inputs$data
  FUN                 <- obj$inputs$fitFUN
  
  this_generation     <- length(obj$generations)
  this_population     <- obj$generations[[this_generation]]$population
  
  n_population        <- length(this_population[,1])
  
  c_fitness           <- mclapply(1:n_population,
                                  ga_fit_itr,
                                  obj,
                                  mc.cores=n_cores)
  
  return(unlist(c_fitness))
}

ga_fit_itr <- function(i,
                       obj,
                       ...){
  
  # use the max_features to modify the ranking fitness of each feature
  n_max_features      <- obj$settings$max_features
  n_fit_penalty       <- obj$settings$fit_penalty
  c_fit_metric        <- obj$settings$fit_metric
  n_generations       <- obj$settings$generations
  v_features          <- obj$inputs$features
  c_predict           <- obj$inputs$predict
  d_data              <- obj$inputs$data
  FUN                 <- obj$inputs$fitFUN
  
  this_generation     <- length(obj$generations)
  this_population     <- obj$generations[[this_generation]]$population
  
  this_gene <- this_population[i,]
  
  this_fea <- v_features[which(this_gene == 1)]
  n_fea <- length(this_fea)
  
  this_fitness <- FUN(this_fea, c_predict, d_data, d_data, ...)
  c_fitness <- this_fitness[[c_fit_metric]] -
    n_fit_penalty * max(0, n_fea - n_max_features)
  
  return(c_fitness)
  
}