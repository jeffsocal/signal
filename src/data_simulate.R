################################################################################
# Copyright (2016) SoCal Bioinformatics Inc. All rights reserved.
# This script is the confidential and proprietary product of SoCal
# Bioinformatics Inc. Any unauthorized reproduction or transfer of
# the contents herein is strictly prohibited.
#
################################################################################
# AUTH:     Jeff Jones | SoCal Bioinofrmatics Inc (SCBI)
# DATE:     2017.02.21
# OWNER:    SoCal Bioinofrmatics Inc
# PROJECT:  SCBI | classifiers
# DESC:     generate simulated data to test classifier algorithms
################################################################################

library(stringr)
library(caret)

data_simulate <- function(
  n_samples=30,
  n_features=5,
  n_sig_features=0,
  data_type='random',
  distribution=rnorm,
  mean_value=1,
  cv_value=0.3,
  diff_value=1.5,
  noise=0.03,
  offset=0.05){
  
  if ( !data_type %in% c('random', 'univariate', 'binary-x', 'binary-positive', 'binary-negative')) {
    stop("data_type must be one of c('random', 'univariate', 'binary-x', 'binary-positive', 'binary-negative')")
  }
  
  # create the matrix
  mat <- matrix(distribution(n_samples*n_features, mean=mean_value, sd=mean_value*cv_value), n_samples, n_features)
  
  # turn into a data frame
  mat <- as.data.frame(mat)
  v_names <- paste0("fea_", str_pad(1:n_features, width=4, pad="0"))
  colnames(mat) <- v_names
  
  # add disease / control status
  v_samples <- 1:n_samples
  v_class_disease <- sample(v_samples, floor(n_samples/2))
  samples_disease <- v_samples[v_class_disease]
  samples_control <- v_samples[-v_class_disease]
  
  mat$patient_status <- 'control'
  mat[v_class_disease,]$patient_status <- 'disease'
  
  v_att <- names(mat)[!grepl('fea', names(mat))]
  
  v_target_features <- c()
  if( data_type == 'univariate' ){
    
    n_sig_features <- max(n_sig_features, 1)
    v_target_features <- sample(v_names, n_sig_features)
    
    for ( i_name in v_target_features ){
      mat[v_class_disease,i_name] <- distribution(length(samples_disease), mean=mean_value, sd=mean_value*cv_value)
      mat[-v_class_disease,i_name] <- distribution(length(samples_control), mean=mean_value*diff_value, sd=mean_value*diff_value*cv_value)
      
    }
    
  } else if( grepl('binary', data_type) ){
    
    n_sig_features <- max(floor(n_sig_features/2)*2, 2)
    v_target_features <- sample(v_names, n_sig_features)
    
    v_noise <- runif(length(mat[,v_target_features[1]]),
                     min=mean_value-mean_value*noise,
                     max=mean_value+mean_value*noise)
    
    mat[,v_target_features[1]] <- distribution(n_samples, mean=mean_value, sd=mean_value*cv_value)
    mat[,v_target_features[2]] <- mat[,v_target_features[1]] + v_noise
    
    # create an offset for parallel distributions
    if( !grepl('x', data_type) ) {
      mat[v_class_disease,v_target_features[2]] <- mat[v_class_disease,v_target_features[2]] + mean_value * offset
      mat[v_class_disease,v_target_features[1]] <- mat[v_class_disease,v_target_features[1]] - mean_value * offset
    }
    
    # pre Process data
    nmat_proc <- preProcess(mat[v_names], method=c('center', 'scale'))
    nmat <- predict(nmat_proc, mat[v_names])
    mat <- cbind(mat[v_att], nmat)
    
    if( grepl('negative', data_type) ) {
      mat[,v_target_features[2]] <- -mat[,v_target_features[2]]
    } else if( grepl('x', data_type) ) {
      mat[v_class_disease,v_target_features[2]] <- -mat[v_class_disease,v_target_features[2]]
    }
    
  }
  
  d_ttest <- ttest_feature(mat, v_features=v_names, c_predict='patient_status', returnTop='sig')
  
  obj_list <- list()
  obj_list[['data']] <- mat
  obj_list[['features']] <- v_names
  obj_list[['target_features']] <- v_target_features
  obj_list[['sig_features']] <- as.character(d_ttest$feature)
  obj_list[['ttest']] <- d_ttest
  
  list_params <- list()
  list_params[['n_samples']] <- n_samples
  list_params[['n_features']] <- n_features
  list_params[['n_sig_features']] <- n_sig_features
  list_params[['data_type']] <- data_type
  #list_params[['distribution']] <- distribution
  list_params[['mean_value']] <- mean_value
  list_params[['cv_value']] <- cv_value
  list_params[['diff_value']] <- diff_value
  list_params[['noise']] <- noise
  list_params[['offset']] <- offset
  
  obj_list[['parameters']] <- list_params
  
  return(obj_list)
}

list_params <- function(obj_list) {
  cat("Fake Data Input Parameters\n")
  l_params <- obj_list[['parameters']]
  for ( i in 1:length(l_params) ){
    cat(" ", str_pad(names(l_params)[i], width=20, side = "right", pad = " "), l_params[[i]], "\n")
  }
}


vector_normalize <- function(v_this, v_other){
  max_other <- max(v_other)
  min_other <- min(v_other)
  
  max_this <- max(v_this)
  min_this <- min(v_this)
  
  v_this <- (v_this - min_this) / (max_this - min_this)
  v_this <- ((v_this) * (max_other - min_other)) + min_other
  return(v_this)
}



features_fromttest <- function(v_features=c(),
                              c_predict='predict',
                              d_data,
                              returnTop=10){
  
  d_fea <- featureTTest(v_features=v_features,
                        c_predict=c_predict,
                        d_data=d_data,
                        returnTop=returnTop)
  
  return(d_fea$feature)
  
}

ttest_feature <- function(d_data,
                         v_features=c(),
                         c_predict='predict',
                         returnTop=10){
  
  predictors <- unique(d_data[,c_predict])
  v_class_a <- which(d_data[,c_predict] == predictors[1])
  v_class_b <- which(d_data[,c_predict] == predictors[2])
  
  v_pvals <- c()
  for ( c_feature in v_features ){
    v_var_a <- d_data[v_class_a, c_feature]
    v_var_b <- d_data[v_class_b, c_feature]
    
    l_ttest <- t.test(v_var_a, v_var_b)
    this_pval <- l_ttest$p.val
    if( is.na(this_pval) | this_pval == 0) {
      this_pval <- 0
    }
    v_pvals <- append(v_pvals, this_pval)
  }
  
  d_pvals <- data.frame(feature=v_features,pval=v_pvals)
  d_pvals <- d_pvals[order(d_pvals$pval),]
  d_pvals$p.adjust.bh <- p.adjust(d_pvals$pval, method='BH')
  
  if( grepl(returnTop,'significant') ){
    return(d_pvals[d_pvals$p.adjust.bh <= 0.05,])
  } else if( is.numeric(returnTop) ) {
    return(d_pvals[1:returnTop,])
  } else {
    stop(paste0("returnTop: '", returnTop, "' not valid, must be 'sig' or as.interger()"))
  }
}

