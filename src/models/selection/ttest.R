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

ttest_sel <- function(v_features=c(),
                      c_predict='predict',
                      d_data=c()){
  
  l_tt <- lapply(1:length(v_features),
                   ttest_wilcox,
                   v_features=v_features,
                   c_predict=c_predict,
                   d_data=d_data)
  
  v_pvals <- unlist(l_tt)
  
  t_pvals <- tibble(
    feature = names(v_pvals),
    p_value = v_pvals,
    p_value_adj = p.adjust(v_pvals, method='BH') 
  ) %>%
    arrange(p_value_adj) %>%
    mutate(rank = row_number())
  
  
  # table of features and values with a rank column  
  return(t_pvals)
  
}

ttest_wilcox <- function(x=1,
                         v_features,
                         c_feature=v_features[x],
                         c_predict,
                         d_data
) {
  
  predictors <- unique(d_data[,c_predict])
  v_class_a <- which(d_data[,c_predict] == predictors[1])
  v_class_b <- which(d_data[,c_predict] == predictors[2])
  
  v_var_a <- d_data[v_class_a, c_feature]
  v_var_b <- d_data[v_class_b, c_feature]
  
  result = tryCatch({
    
    l_ttest <- wilcox.test(v_var_a, v_var_b)
    this_pval <- l_ttest$p.val
    
    if( is.na(this_pval) | this_pval == 1) {
      this_pval <- 1
    }
    
    l <- list()
    l[[c_feature]] <- this_pval
    return(l)
  }, warning = function(w) {
    l <- list()
    l[[c_feature]] <- 1
    return(l)
  }, error = function(e) {
    l <- list()
    l[[c_feature]] <- 1
    return(l)
  }, finally = {
  })
  
  
}

# ttest_wilcox <- function(x=1,
#                          v_features,
#                          c_feature=v_features[x],
#                          c_predict,
#                          d_data
# ) {
#   
#   predictors <- unique(d_data[,c_predict])
#   v_class_a <- which(d_data[,c_predict] == predictors[1])
#   v_class_b <- which(d_data[,c_predict] == predictors[2])
#   
#   v_var_a <- d_data[v_class_a, c_feature]
#   v_var_b <- d_data[v_class_b, c_feature]
#   
#   result = tryCatch({
#     l_ttest <- wilcox.test(v_var_a, v_var_b)
#     this_pval <- l_ttest$p.val
#     
#     if( is.na(this_pval) | this_pval == 1) {
#       this_pval <- 1
#     }
#     
#     names(this_pval) <- c_feature
#     return(this_pval)
#   }, warning = function(w) {
#     this_pval <- 1
#     names(this_pval) <- c_feature
#     return(this_pval)
#   }, error = function(e) {
#     this_pval <- 1
#     names(this_pval) <- c_feature
#     return(this_pval)
#   }, finally = {
#     this_pval <- 1
#     names(this_pval) <- c_feature
#     return(this_pval)
#   })
#   
#   
# }
