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

# INPUT  : object.1 LIST of MATRIX.populations and VECTOR.fitness
#        : object.2 VECTOR of MATRIX indices (not sure why, but could be useful)
#        : object.3 FUNCTION
#
# OUTPUT : object LIST of MATRIX.populations and VECTOR.fitness

# public
ga_selection <- function(obj, 
                        parents=NULL, 
                        FUN=selSimpleRank, ...){
  
  # use the max_features to modify the ranking fitness of each feature
  n_max_features      <- obj$settings$max_features
  n_generations       <- obj$settings$generations
  v_features          <- obj$inputs$features
  c_predict           <- obj$inputs$predict
  d_data              <- obj$inputs$data
  
  this_generation     <- length(obj$generations)-1
  this_population     <- obj$generations[[this_generation]]$population
  this_fitness        <- obj$generations[[this_generation]]$fitness
  
  n_population        <- length(this_population[,1])
  n_max_population    <- min(n_population, obj$settings$max_population)
  
  sel <- FUN(-this_fitness, n_max_population)
  
  out <- list(population = this_population[sel,,drop=FALSE],
              fitness = this_fitness[sel])
  
  return(out)
}




# SELECTION

# private
probSelection <- function(p_fitness,
                          max_population=10){
  
  n_population <- length(p_fitness)
  
  v_selection <- sample(1:n_population, size = max_population, 
                        prob = pmin(pmax(0, p_fitness), 1, na.rm = TRUE),
                        replace = TRUE)
  return(v_selection)
}

# private
selSimpleRank <- function(v_fitness,
                             max_population=10){
  
  v_selection <- order(v_fitness)[1:max_population]
  return(v_selection)
}

# private
selNonlinearRank <- function(v_fitness,
                             max_population=10){
  
  # Michalewicz (1996) Genetic Algorithms + Data Structures = Evolution Programs. p. 60
  r_values <- rank(v_fitness, ties.method = 'random')
  q <- 2/length(r_values)
  r <- 2/(length(r_values) * (length(r_values)-1))
  p_values <- q - (r_values-1) * r
  
  v_selection <- probSelection(p_values, max_population)
  return(v_selection)
}


# private
selExponentialRank <- function(v_fitness,
                               max_population=10){
  # Michalewicz (1996) Genetic Algorithms + Data Structures = Evolution Programs. p. 60
  r_values <- rank(v_fitness, ties.method = 'random')
  q <- 0.25
  p_values <- q*(1-q)^(r_values-1)
  
  v_selection <- probSelection(p_values, max_population)
  return(v_selection)
}


# private
selFitnessProp <- function(v_fitness,
                           max_population=10){
  # roulette wheel selection (fitness proportionate selection)
  # http://en.wikipedia.org/wiki/Fitness_proportionate_selection
  v_fitness <- v_fitness - min(v_fitness, na.rm=T)
  p_values <- v_fitness/sum(v_fitness)
  
  v_selection <- probSelection(p_values, max_population)
  return(v_selection)
}


# private
selTournament <- function(v_fitness, 
                          max_population=10,
                          k = 3) {
  
  # (unbiased) Tournament selection 
  v_selection <- rep(NA, max_population)
  
  for(i in 1:max_population){ 
    s <- sample(1:max_population, size = k)
    v_selection[i] <- s[which.max(v_fitness[s])]
  }
  
  return(v_selection)
}

