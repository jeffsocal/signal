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

# requires
library(parallel)
library(caret)


# data obj model
#
# obj                              list
# .......inputs                        list
# ...............dData                  matrix / data.frame
# ...............cPredictor             character (column in dData)
# ...............vFeatures              vector
# ...............mSeed                  matrix / data.frame
# ...............fFitness               function
# .......settings                      list
# ...............nGenes                 int (features)
# ...............nPopulation            int (feature sets)
# ...............nGenerations           int
# ...............nMaxActiveGenes        int (target n features)
# ...............nFitPenalty            int
# ...............nFitMetric             int
# ...............pSurvive               float (probability dist of survivors)
# ...............pCrossover             float
# ...............pMutation              float
# ...............pSwap                  float
# ...............pShuffle               float
# .......generations                   list
# ..................[[i]]               list
# .......................population      matrix
# .......................fitness         vector
# .......solution                      list
# ...............features               vector
# ...............fitness                numeric



# public
ga <- function(obj,
               verbose=T,
               ...) {
  
  if( verbose == T )  cat("\n")
  
  n_generations       <- obj$settings$generations
  m_seed              <- obj$inputs$seed
  
  for ( gen in 1:n_generations ){
    
    if( verbose == T )
      cat(str_pad(paste("   Generation", gen, " "), width=20))
    
    obj$generations[[gen]] <- list()    
    
    if( gen == 1 ) {
      obj$generations[[gen]]$population <- m_seed
      obj$generations[[gen]]$fitness <- ga_fit(obj)
      if( verbose == T )
        ga_gen_stats(obj$generations[[gen]]$fitness)
      
      next()
    } 
    
    last_gen <- new_gen <- ga_selection(obj)
    n_last_gen_size <- length(last_gen$population[,1])
    # --- do stuff to the population
    for ( gac in 2:min(10, n_last_gen_size )){
      new_gen$population <- rbind(new_gen$population,
                                  gaCrossover(last_gen, 1:gac))
    }
    
    new_gen$population <- rbind(new_gen$population,
                                gaMutate(last_gen))
    
    new_gen$population <- rbind(new_gen$population,
                                gaSwap(last_gen))
    
    new_gen$population <- rbind(new_gen$population,
                                gaShuffle(last_gen))
    
    obj$generations[[gen]] <- new_gen
    obj$generations[[gen]] <- gaDedup(new_gen)
    
    obj$generations[[gen]]$fitness <- ga_fit(obj)
    
    if( verbose == T )
      ga_gen_stats(obj$generations[[gen]]$fitness)
    
    # EARLY EXIT ====================================================
    # immediate last gen shows little change
    this_auc <- ga_best(obj,gen,verbose=F)
    last_auc <- ga_best(obj,gen-1,verbose=F)
    if( abs(this_auc$max - last_auc$max) < 0.0001 &
        (abs(this_auc$max - last_auc$mean) / this_auc$max) < 0.025 )
      return(obj)
    # last 10 gens show no change in max AUC
    if( length(obj$generations) > 10 ) {
      last_auc <- ga_best(obj,gen-10,verbose=F)  
      if( abs(this_auc$max - last_auc$max) < 0.0001 )
        return(obj)
    }
    # EARLY EXIT ====================================================
    
  }
  
  return(obj)
}

ga_gen_stats <- function(this_fit){
  cat(str_pad(paste("n", length(this_fit)), width=7),
      str_pad(paste("max", round(max(this_fit, na.rm=T), 3)), width=12),
      str_pad(paste("min", round(min(this_fit, na.rm=T), 3)), width=12),
      str_pad(paste("mean", round(mean(this_fit, na.rm=T), 3)), width=12),
      "\n")
}



# public
ga_initalize <- function(d_data,
                         c_predict,
                         v_features, 
                         m_seed=c(),
                         n_cores=1,
                         n_max_generations=25,
                         n_max_population=50,
                         n_max_features=ceiling(length(v_features)*.05),
                         n_fit_penalty=0.01,
                         c_fit_metric='auc',
                         f_fitness=svm_cls,
                         ...
) {
  
  obj <- list()
  
  # .......inputs                                   list
  obj[['inputs']] <- list(
    features = v_features,
    predict = c_predict,
    data = d_data,
    seed = m_seed,
    fitFUN = f_fitness
  )
  
  # .......settings                                 list
  obj[['settings']] <- list(
    cores = n_cores,
    genes = length(v_features),
    generations = n_max_generations,
    max_population = n_max_population,
    max_features = n_max_features,
    fit_penalty = n_fit_penalty,
    fit_metric = c_fit_metric,
    prob_survive = 1,
    prob_crossover = 1,
    prob_mutation = 1,
    prob_swap = 1,
    prob_shuffle = 1
  )
  
  # if( length(m_seed) == 0) {
  #   m_seed <- ga_seed(obj, method='rand', n_comb=2, limit=5e5)
  #   obj$inputs$seed <- m_seed
  # }
  
  # .......generations                   list
  obj[['generations']] <- list()
  # ..................[[i]]               list
  # .......................population      matrix
  # .......................fitness         vector
  
  # .......solution                      list
  obj[['solution']] <- list()
  # ...............features               vector
  # ...............fitness                numeric
  
  return(obj)
}


ga_seed <- function(obj, 
                    method='bf',
                    n_comb=2,
                    limit=NULL,
                    ...) {
  
  if( is.null(limit) ) 
    limit <- 100
  
  v_features          <- obj$inputs$features
  n_features          <- length(v_features)
  
  n_population        <- choose(n_comb:n_features, n_comb)
  n_population        <- max(which(n_population <= limit))
  m_use               <- t(combn(sample(1:n_features, n_population), n_comb))
  
  m_seed <- matrix(0, ncol=n_features, nrow=dim(m_use)[1])
  
  for( i in 1:n_population ) {
    m_seed[i,m_use[i,1]] <- 1
    m_seed[i,m_use[i,2]] <- 1
  }
  
  return(m_seed)
}



ga_update <- function(obj, ...) {
  
  latest_gen <- obj$generations[[length(obj$generations)]]
  i_fittest  <- which(latest_gen$fitness == max(latest_gen$fitness, na.rm=T))[1]
  
  obj$solution[['features']] <- latest_gen$population[i_prog_fit,]
  obj$solution[['fitness']]  <- latest_gen$fitness[i_prog_fit]
  
}


ga_summary <- function(obj, ...) {
  
  d_settings <- as.data.frame(unlist(obj$settings))
  d_settings$variable <- rownames(d_settings)
  colnames(d_settings)[1] <- 'input'
  
  cat("\nGA Feature Select Input Parameters\n")
  print(d_settings[,c('variable','input')], row.names=F)
  
}