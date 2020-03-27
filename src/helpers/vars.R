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

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  as.character(tolower(os))
}

# SET VARIABLES
################################################################################

data_path <- "Z:/data"
if(get_os() == 'linux'){
  data_path <- "~/qnapd/data"
}

################################################################################