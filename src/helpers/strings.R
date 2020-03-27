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

library(stringr)

all_messages <<- c()

print_message <- function(str_message='',str_value='NULL'){
  
  #Sys.sleep(1)
  date_stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  
  if(str_value == 'NULL'){
    str_message <- paste(str_message, "\n")
  } else {
    str_message <- paste(str_pad(str_message, 50, side = "right", pad = " "), ":", str_value, "\n", sep="")
  }
  cat(str_message)
  
  all_messages[[length(all_messages)+1]] <<- str_message
}

print_div <- function(){
  cat("\n============================================================\n")
}