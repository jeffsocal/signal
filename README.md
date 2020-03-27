# signal
A classification R project

This package utilizes various feature discovery and classification schemes

### Feature Discovery 

The top_n fetaures will be passed forward to classification.

#### t.test
A Wilcox ttest based on class designation.

#### enet
Using ElasticNet package R::elasticnet.

#### rf
Using RandomForest package R::randomForest.

#### mine
Using Maximal Information Non-parametric Exploration package R::minerva.

#### bf
Using a brute-force all possible combination SVM auc ranked.


### Classification

#### svm, ksvm
An SVM classification algorithm based on the R packages.

 - R::e1071 for svm 6 combinations
 -- types: eps-regression, nu-regression
 -- kernels: linear, radial, polynomial
 - R::kernlab for svm 81 combinations
 -- types: C-svc, nu-svc, C-bsvc, spoc-svc, kbb-svc, one-svc, eps-svr, nu-svr, eps-bsvr
 -- kernels: rbfdot, polydot, vanilladot, tanhdot, laplacedot, besseldot, anovadot, splinedot, stringdot

#### glm
Using generalized linear model.

#### rf
Using RandomForest package R::randomForest.


## Setup and exection

#### User inputs
- modeling_data: the data matrix 
- usr_project: name the project
- usr_dataset: name the data set
- usr_n_features: n-features to target, string: '1', '1:20', '2,5,10'
- usr_features: down select the feaures as vector or designate 'all'
- usr_cfv_split: define how the folds are split, null = random, by column name for paring
- usr_method: define the selection method
- usr_model: define the classification model
- usr_n_fold: number of folds
- usr_n_reps: number of replicates
- usr_n_cores: parallel processing cores
- usr_verbose: print status messages for logging


##### example
`
# import the code
source("./src/signal.R")

# define the data
r_data <- readRDS("./dat/000_unittest/univariate_100fea-100samples_2sig.R")

# create an integer based class designation
r_data$data[,"patient_integer"] <- as.numeric(r_data$data[,'patient_status'] == "disease")

# apply preprocessing 
modeling_data <- preprocess(r_data$data, fea, "patient_integer")

class_model <- signal(
  modeling_data,
  usr_project = 'unittest',
  usr_dataset = 'unittest_100x100',
  usr_n_features = '2',
  usr_features = 'all',
  usr_cfv_split = null,
  usr_method = 'ttest',
  usr_model = 'svm05',
  usr_n_fold = 10,
  usr_n_reps = 10,
  usr_n_cores = 1,
  usr_verbose = T
)
`