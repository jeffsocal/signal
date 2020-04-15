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

source("./src/models/genealg/gafeat.R")
source("./src/models/genealg/ga_explore.R")
source("./src/models/genealg/ga_fit.R")
source("./src/models/genealg/ga_manipulations.R")
source("./src/models/genealg/ga_selections.R")

ga_sel <- function(v_features,
                   c_predict,
                   d_data,
                   n_features=10,
                   n_cores=1,
                   n_max_generations=20,
                   n_fit_penalty=0.02,
                   seed_method='rand',
                   seed_n_comb=3,
                   seed_n_limit=1e3,
                   verbose=T
) {
  
  obj <- ga_initalize(v_features=v_features,
                      c_predict=c_predict,
                      d_data=d_data,
                      m_seed=c(),
                      n_cores=n_cores,
                      n_max_generations=n_max_generations,
                      n_max_features=n_features,
                      n_fit_penalty=n_fit_penalty
  )
  
  # Create a seed after the object initalization
  obj$inputs$seed <- ga_seed(obj,
                             method=seed_method,
                             n_comb=seed_n_comb,
                             limit=seed_n_limit)
  
  
  obj <- ga(obj,
            verbose=verbose)
  
  last_gen <- obj$generations[[length(obj$generations)]]$population %>% 
    as_tibble()
  
  colnames(last_gen) <- obj$inputs$features
  
  t_out <- last_gen %>% 
    gather(., key="feature", value="frequency") %>% 
    group_by(feature) %>% 
    summarise(frequency = sum(frequency)) %>% 
    arrange(desc(frequency)) %>% 
    mutate(rank = row_number())
  
  return(t_out)
  
}

