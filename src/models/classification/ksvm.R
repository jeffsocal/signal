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

library(e1071)
library(ROCR)
library(kernlab)

ksvmt <- function(v_features,
                  c_predict,
                  d_train,
                  d_test,
                  verbose=F,
                  type='eps-svr',
                  kernel='anovadot'){

  v_features <- v_features[!is.na(v_features)]

  out <- list(pred=NULL,
              auc=0,
              fpr=NULL,
              tpr=NULL,
              features=v_features,
              predictor=c_predict,
              prediction=NULL,
              observation=NULL)

  if(length(v_features) == 0)
    return(out)

  if( grepl('svc', type ) ) {
    d_train[,c_predict] <- as.factor(d_train[,c_predict])
    d_test[,c_predict] <- as.factor(d_test[,c_predict])
    # removed to insead use a predefined numeric class class
    #   } else {
    #     d_train[,c_predict] <- as.numeric(d_train[,c_predict] == 'disease')
    #     d_test[,c_predict] <- as.numeric(d_test[,c_predict] == 'disease')
  }


  out_model <- exe_ksvm(d_tr=d_train,
                        d_ts=d_test,
                        c_pd=c_predict,
                        v_fe=v_features,
                        type=type,
                        kernel=kernel)

  if( is.null(out_model$prediction) )
    return(out)

  #calculate the ROC AUC value
  pred     				<- prediction(out_model$prediction, d_test[,c_predict])
  perf		 				<- performance(pred, "auc")

  out[['pred']] <- pred
  out[['perf']] <- perf
  out[['auc']] <- unlist(perf@y.values)

  perf    	 				<- performance(pred, "tpr", "fpr")
  out[['fpr']] <- unlist(perf@x.values)
  out[['tpr']] <- unlist(perf@y.values)

  out[['classifier']] <- out_model$classifier
  out[['prediction']] <- out_model$prediction
  out[['observation']] <- d_test[,c_predict]

  if( verbose == T) {
    cat("--> KSVM:", round(out[['auc']],3), ":")
    str(v_features)
  }

  return(out)
}


exe_ksvm <- function(d_tr,
                     d_ts,
                     c_pd,
                     v_fe,
                     type,
                     kernel){

  result = tryCatch({

    # formula
    flma <- as.formula(paste0(c_pd, "~."))
    capture.output(svm_model <- ksvm(flma, data=d_tr[,c(v_fe, c_pd)], type=type, kernel=kernel, na.action = na.fail, scaled = F))

    # matrix
    # d_tr <- as.matrix(d_tr)
    # d_ts <- as.matrix(d_ts)
    # svm_model <- ksvm(d_tr[,v_fe], d_tr[,c_pd], type=type, kernel=kernel)

    predict_model <- as.vector(predict(svm_model, d_ts[,c(v_fe)]))
    # removed to insead use a predefined numeric class class
    #     if( grepl('svc', type )) {
    #       predict_model <- as.numeric(predict_model == 'disease')
    #     }
    names(predict_model) <- rownames(d_ts)

    out_model <- list()
    out_model$classifier <- svm_model
    out_model$prediction <- predict_model
    #calculate the ROC AUC value
    return(out_model)

  }, warning = function(w) {
    #stop(w)
    return(NULL)
  }, error = function(e) {
    #stop(e)
    return(NULL)
  }, finally = {

  })
}



ksvm_models <- function(){

  #   type <- c('C-svc','nu-svc','C-bsvc','spoc-svc','kbb-svc','one-svc','eps-svr','nu-svr','eps-bsvr')
  #   kernel <- c('rbfdot','polydot','vanilladot','tanhdot','laplacedot','besseldot','anovadot','splinedot','stringdot')
  #
  #   combs <- expand.grid(type, kernel)
  #   for ( i in 1:dim(combs)[1] ){
  #     cat(paste0("l_models[['ksvm", str_pad(i, width=2, pad=0), "']] <- function(...) ksvmt(..., type='", combs[i,1], "', kernel='", combs[i,2], "')"),"\n")
  #   }

  l_models <- list()
  l_models[['ksvm01']] <- function(...) ksvmt(..., type='C-svc', kernel='rbfdot')
  l_models[['ksvm02']] <- function(...) ksvmt(..., type='nu-svc', kernel='rbfdot')
  l_models[['ksvm03']] <- function(...) ksvmt(..., type='C-bsvc', kernel='rbfdot')
  l_models[['ksvm04']] <- function(...) ksvmt(..., type='spoc-svc', kernel='rbfdot')
  l_models[['ksvm05']] <- function(...) ksvmt(..., type='kbb-svc', kernel='rbfdot')
  l_models[['ksvm06']] <- function(...) ksvmt(..., type='one-svc', kernel='rbfdot')
  l_models[['ksvm07']] <- function(...) ksvmt(..., type='eps-svr', kernel='rbfdot')
  l_models[['ksvm08']] <- function(...) ksvmt(..., type='nu-svr', kernel='rbfdot')
  l_models[['ksvm09']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='rbfdot')
  l_models[['ksvm10']] <- function(...) ksvmt(..., type='C-svc', kernel='polydot')
  l_models[['ksvm11']] <- function(...) ksvmt(..., type='nu-svc', kernel='polydot')
  l_models[['ksvm12']] <- function(...) ksvmt(..., type='C-bsvc', kernel='polydot')
  l_models[['ksvm13']] <- function(...) ksvmt(..., type='spoc-svc', kernel='polydot')
  l_models[['ksvm14']] <- function(...) ksvmt(..., type='kbb-svc', kernel='polydot')
  l_models[['ksvm15']] <- function(...) ksvmt(..., type='one-svc', kernel='polydot')
  l_models[['ksvm16']] <- function(...) ksvmt(..., type='eps-svr', kernel='polydot')
  l_models[['ksvm17']] <- function(...) ksvmt(..., type='nu-svr', kernel='polydot')
  l_models[['ksvm18']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='polydot')
  l_models[['ksvm19']] <- function(...) ksvmt(..., type='C-svc', kernel='vanilladot')
  l_models[['ksvm20']] <- function(...) ksvmt(..., type='nu-svc', kernel='vanilladot')
  l_models[['ksvm21']] <- function(...) ksvmt(..., type='C-bsvc', kernel='vanilladot')
  l_models[['ksvm22']] <- function(...) ksvmt(..., type='spoc-svc', kernel='vanilladot')
  l_models[['ksvm23']] <- function(...) ksvmt(..., type='kbb-svc', kernel='vanilladot')
  l_models[['ksvm24']] <- function(...) ksvmt(..., type='one-svc', kernel='vanilladot')
  l_models[['ksvm25']] <- function(...) ksvmt(..., type='eps-svr', kernel='vanilladot')
  l_models[['ksvm26']] <- function(...) ksvmt(..., type='nu-svr', kernel='vanilladot')
  l_models[['ksvm27']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='vanilladot')
  l_models[['ksvm28']] <- function(...) ksvmt(..., type='C-svc', kernel='tanhdot')
  l_models[['ksvm29']] <- function(...) ksvmt(..., type='nu-svc', kernel='tanhdot')
  l_models[['ksvm30']] <- function(...) ksvmt(..., type='C-bsvc', kernel='tanhdot')
  l_models[['ksvm31']] <- function(...) ksvmt(..., type='spoc-svc', kernel='tanhdot')
  l_models[['ksvm32']] <- function(...) ksvmt(..., type='kbb-svc', kernel='tanhdot')
  l_models[['ksvm33']] <- function(...) ksvmt(..., type='one-svc', kernel='tanhdot')
  l_models[['ksvm34']] <- function(...) ksvmt(..., type='eps-svr', kernel='tanhdot')
  l_models[['ksvm35']] <- function(...) ksvmt(..., type='nu-svr', kernel='tanhdot')
  l_models[['ksvm36']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='tanhdot')
  l_models[['ksvm37']] <- function(...) ksvmt(..., type='C-svc', kernel='laplacedot')
  l_models[['ksvm38']] <- function(...) ksvmt(..., type='nu-svc', kernel='laplacedot')
  l_models[['ksvm39']] <- function(...) ksvmt(..., type='C-bsvc', kernel='laplacedot')
  l_models[['ksvm40']] <- function(...) ksvmt(..., type='spoc-svc', kernel='laplacedot')
  l_models[['ksvm41']] <- function(...) ksvmt(..., type='kbb-svc', kernel='laplacedot')
  l_models[['ksvm42']] <- function(...) ksvmt(..., type='one-svc', kernel='laplacedot')
  l_models[['ksvm43']] <- function(...) ksvmt(..., type='eps-svr', kernel='laplacedot')
  l_models[['ksvm44']] <- function(...) ksvmt(..., type='nu-svr', kernel='laplacedot')
  l_models[['ksvm45']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='laplacedot')
  l_models[['ksvm46']] <- function(...) ksvmt(..., type='C-svc', kernel='besseldot')
  l_models[['ksvm47']] <- function(...) ksvmt(..., type='nu-svc', kernel='besseldot')
  l_models[['ksvm48']] <- function(...) ksvmt(..., type='C-bsvc', kernel='besseldot')
  l_models[['ksvm49']] <- function(...) ksvmt(..., type='spoc-svc', kernel='besseldot')
  l_models[['ksvm50']] <- function(...) ksvmt(..., type='kbb-svc', kernel='besseldot')
  l_models[['ksvm51']] <- function(...) ksvmt(..., type='one-svc', kernel='besseldot')
  l_models[['ksvm52']] <- function(...) ksvmt(..., type='eps-svr', kernel='besseldot')
  l_models[['ksvm53']] <- function(...) ksvmt(..., type='nu-svr', kernel='besseldot')
  l_models[['ksvm54']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='besseldot')
  l_models[['ksvm55']] <- function(...) ksvmt(..., type='C-svc', kernel='anovadot')
  l_models[['ksvm56']] <- function(...) ksvmt(..., type='nu-svc', kernel='anovadot')
  l_models[['ksvm57']] <- function(...) ksvmt(..., type='C-bsvc', kernel='anovadot')
  l_models[['ksvm58']] <- function(...) ksvmt(..., type='spoc-svc', kernel='anovadot')
  l_models[['ksvm59']] <- function(...) ksvmt(..., type='kbb-svc', kernel='anovadot')
  l_models[['ksvm60']] <- function(...) ksvmt(..., type='one-svc', kernel='anovadot')
  l_models[['ksvm61']] <- function(...) ksvmt(..., type='eps-svr', kernel='anovadot') #<-- best regression
  l_models[['ksvm62']] <- function(...) ksvmt(..., type='nu-svr', kernel='anovadot')
  l_models[['ksvm63']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='anovadot')
  l_models[['ksvm64']] <- function(...) ksvmt(..., type='C-svc', kernel='splinedot')
  l_models[['ksvm65']] <- function(...) ksvmt(..., type='nu-svc', kernel='splinedot')
  l_models[['ksvm66']] <- function(...) ksvmt(..., type='C-bsvc', kernel='splinedot')
  l_models[['ksvm67']] <- function(...) ksvmt(..., type='spoc-svc', kernel='splinedot')
  l_models[['ksvm68']] <- function(...) ksvmt(..., type='kbb-svc', kernel='splinedot')
  l_models[['ksvm69']] <- function(...) ksvmt(..., type='one-svc', kernel='splinedot')
  l_models[['ksvm70']] <- function(...) ksvmt(..., type='eps-svr', kernel='splinedot')
  l_models[['ksvm71']] <- function(...) ksvmt(..., type='nu-svr', kernel='splinedot')
  #   l_models[['ksvm72']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='splinedot')
  #   l_models[['ksvm73']] <- function(...) ksvmt(..., type='C-svc', kernel='stringdot')
  #   l_models[['ksvm74']] <- function(...) ksvmt(..., type='nu-svc', kernel='stringdot')
  #   l_models[['ksvm75']] <- function(...) ksvmt(..., type='C-bsvc', kernel='stringdot')
  #   l_models[['ksvm76']] <- function(...) ksvmt(..., type='spoc-svc', kernel='stringdot')
  #   l_models[['ksvm77']] <- function(...) ksvmt(..., type='kbb-svc', kernel='stringdot')
  #   l_models[['ksvm78']] <- function(...) ksvmt(..., type='one-svc', kernel='stringdot')
  l_models[['ksvm79']] <- function(...) ksvmt(..., type='eps-svr', kernel='stringdot')
  l_models[['ksvm80']] <- function(...) ksvmt(..., type='nu-svr', kernel='stringdot')
  #   l_models[['ksvm81']] <- function(...) ksvmt(..., type='eps-bsvr', kernel='stringdot')


  return(l_models)
}
