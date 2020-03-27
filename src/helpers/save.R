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


saveResultObject <- function(object,
                             project,
                             dataset=".",
                             process='generic'){

  date_stamp <- format(Sys.time(), "%Y%m%d%H%M%S")

  path <- "./dat/_results/"
  dir.create(file.path(path, project))
  path <- paste0(path, project, "/")

  saveRDS(object, file=paste0(path,
                              date_stamp, "_",
                              project, "_",
                              dataset, "_",
                              process, ".rds"))
}

