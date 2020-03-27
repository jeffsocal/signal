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

# set up local environment and load code scripts
rm(list=ls(all=TRUE))
source("./src/data_simulate.R")

#===============================================================================

#===============================================================================
# Create completely random normaly-distributed data with rnorm
f_data <- data_simulate(data_type='random',
                         n_features = 50,
                         n_samples = 50
                         )

# pop out the data.frame
data_frame <- f_data$data

# pop out a vector of features
v_features <- f_data$features

# univariate t.test has alreday been done, take a look
# only BH adjusted significant features are listed
print(f_data$ttest)

# pop out a vector of significant features, if they exist
v_significant_features  <- f_data$sig_features

# list the input parameters
list_params(f_data)



#===============================================================================
# Create a random normaly-distributed data set with a single significant feature
f_data <- data_simulate(data_type='univariate',
                         n_features = 100,
                         n_samples = 100,
                         n_sig_features=2
                         )

saveRDS(f_data, "./dat/000_unittest/univariate_100fea-100samples_2sig.R")



#===============================================================================
# Create a random normaly-distributed data set with 2 interacting
# features, with positive association
f_data <- data_simulate(data_type   ='binary-x',
                         n_features  = 5000,
                         n_samples   = 250,
                         mean_value  = 10,
                         cv_value    = 0.2,
                         diff_value  = 1.5,
                         noise       = 0.2,
                         offset      = 0.0)

# works only for two features, in this case the target_features
# plot the two target_features (i.e. the features chosen to interact)

dat <- f_data$data
fea <- f_data$target_features
ggplotBinaryInteraction(dat, fea)
ggplotVolcanoPlot(dat, fea)



f_data <- data_simulate(data_type='univariate',
                         n_features = 2,
                         n_samples = 300)

# works only for two features, in this case the target_features
# plot the two target_features (i.e. the features chosen to interact)
dat <- f_data$data
dat <- dat[order(dat$patient_status),]
fea <- f_data$target_features
ggplotBinaryInteraction(dat, c(fea,f_data$features[2]))

#===============================================================================
# Create a random normaly-distributed data set with 2 interacting
# features, with negative association
f_data <- data_simulate(data_type='binary-negative',
                         n_features = 50,
                         n_samples = 50)

# plot
dat <- f_data$data
fea <- f_data$target_features
ggplotBinaryInteraction(dat, fea)




#===============================================================================
# Create a random normaly-distributed data set with 2 interacting
# features, with disease / control alternate associations
f_data <- data_simulate(data_type='binary-x',
                         n_features = 50,
                         n_samples = 100)

# plot
dat <- f_data$data
fea <- f_data$target_features
ggplotBinaryInteraction(dat, fea)





