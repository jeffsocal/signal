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

library(pryr)

memsize <- function(dig=2){
  
  x <- as.numeric(mem_used()) / 10^9
  
  return(paste(round(x,dig), "Gb"))
}
