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
#        : object.2 VECTOR of MATRIX indices
#
# OUTPUT : object MATRIX.populations
#        : (not including the INPUT)

# public
gaDedup <- function(object, 
                    parents=NULL, 
                    ...) {
  
  if( is.null(parents) )
    parents <- 1:length(object$population[,1])
  
  m_prnt <- m_chld <- as.matrix(object$population[parents,])
  n_col <- ncol(m_prnt)
  n_row <- nrow(m_prnt)
 
  v_chr <- as.character(
            unlist(
              lapply(
                as.data.frame(t(m_prnt)), 
                paste, collapse="")
              )
            )
  
  v_zro <- as.numeric(
            unlist(
              lapply(
                as.data.frame(t(m_prnt)), sum)
              )
            )
  
  v_dup <- union(
            which(duplicated(v_chr)), 
            which(v_zro == 0)
            )
  
  if(length(v_dup) > 0){
    object$population <- object$population[-v_dup,]
    object$fitness <- object$fitness[-v_dup]
  }
  
  return(object)
}

# public
gaMutate <- function(object, 
                     parents=NULL, 
                     method=runif,
                     ...) {
  
  if( is.null(parents) )
    parents <- 1:length(object$population[,1])
  
  m_prnt <- m_chld <- as.matrix(object$population[parents,])
  n_col <- ncol(m_prnt)
  n_row <- nrow(m_prnt)
  
  # TODO:
  # take into consideration the population distribution and
  # frequency of observation for each gene, mutating IN a
  # gene that otherwise was never selected 
  
  # uniform random mutation
  for ( i in 1:n_row ){    
    j <- sample(1:n_col, size = 1)
    m_chld[i,j] <- abs(m_chld[i,j]-1)
  }
  
  return(m_chld)
}


# public
gaSwap <- function(object, 
                   parents=NULL, 
                   ...) {
  
  if( is.null(parents) )
    parents <- 1:length(object$population[,1])
  
  # a mechanism by which single alleles are swapped
  # i.e. for a [1] -> 0 also perform a [0] -> 1
  
  m_prnt <- m_chld <- object$population[parents,]
  n_col <- ncol(m_prnt)
  n_row <- nrow(m_prnt)

  # uniform random swap
  for ( i in 1:n_row ){    
    for ( k in 0:1 ){
      w <- which(m_prnt[i,] == k)
      if( length(w) == 0 ) break()
      j <- sample(w, size = 1)
      m_chld[i,j] <- abs(m_chld[i,j]-1)
    }
  }
  
  return(m_chld)
}


# public
gaShuffle <- function(object, 
                      parents=NULL, 
                      ...) {
  
  if( is.null(parents) )
    parents <- 1:length(object$population[,1])
  
  # uniform random shuffle
  
  m_prnt <- m_chld <- as.matrix(object$population[parents,])
  n_col <- ncol(m_prnt)
  n_row <- nrow(m_prnt)
  
  v_val <- as.vector(m_prnt)
  n_val <- length(m_prnt)
  m_chld <- matrix(sample(v_val, n_val, replace=F), ncol=n_col, nrow=n_row)
  
  return(m_chld)
}


# public
gaCrossover <- function(object, 
                        parents=NULL, 
                        ...) {
  
  if( is.null(parents) )
    parents <- 1:length(object$population[,1])
  
  # Single-point crossover
  # ensure there are only 2 parents
  parents <- sample(parents, size = 2, replace = T)
  
  m_prnt <- m_chld <- as.matrix(object$population[parents,])
  n_col <- ncol(m_prnt)
  n_row <- nrow(m_prnt)
  
  # ensure there are >= 2 genes
  if ( n_row <= 1 )
    return(m_prnt)
  
  m_chld <- matrix(as.double(NA), nrow = n_row, ncol = n_col)
  
  i_crossover <- sample(1:(n_col-1), size = 1)

  v_ls <- 1:i_crossover       # left side
  v_rs <- (i_crossover+1):n_col   # right side
  
  m_chld[1,] <- c( m_prnt[1,v_ls], m_prnt[2,v_rs] )
  m_chld[2,] <- c( m_prnt[2,v_ls], m_prnt[1,v_rs] )
  
  return(m_chld)
}

