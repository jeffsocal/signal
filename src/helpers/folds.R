################################################################################
# Copyright (2016) SoCal Bioinformatics Inc. All rights reserved.
# This script is the confidential and proprietary product of SoCal 
# Bioinformatics Inc. Any unauthorized reproduction or transfer of 
# the contents herein is strictly prohibited.
#
################################################################################
# AUTH:     Jeff Jones | SoCal Bioinofrmatics Inc
# DATE:     2016.11.16
# OWNER:    SoCal Bioinofrmatics Inc
################################################################################

# mimic caret
createFoldsFromID <- function(x, k = 10, list = TRUE, returnTrain = FALSE){
  ux <- sample(unique(x))
  lx <- length(ux)
  sx <- split(ux, sort(1:lx%%k))
  ox <- list()
  for( i in 1:k ){
    wx <- which( x %in% sx[[i]])
    ox[[paste0("Fold", str_pad(i, 2, 'left', 0))]] <- wx 
  }
  return(ox)
}
