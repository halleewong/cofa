## --- Summary --------------------------------------------------------------
# Functions for training and testing (c-statistic) ensemble models

## --- Load packages --------------------------------------------------------

library(fifer)
library(glmnet)
library(caret)

## --- makeGroups --------------------------------------------------------

makeGroups <- function(df, colName, threshold){
  # parameters:
  #   df - a data.frame containing colName
  #   colName - the column to group 
  #   threshold - all levels of colName less than this threshold will be consolidates
  #               into a "other" category
  # returns:
  #   a vector with the new levels of the categorical variable for each case
  
  # get names of chapter that are too small
  tbl = data.frame(table(df[,colName]))
  tooSmall = tbl[tbl$Freq < threshold,"Var1"]
  
  # convert to factor
  newCol = as.factor(df[,colName])
  
  # bag all small chapters into "other"
  levels(newCol)[levels(newCol) %in% tooSmall] <- "other"
  newCol[is.na(newCol)] <- factor("other")
  
  return(newCol)
} 

## --- makeTrainTest --------------------------------------------------------

makeTrainTest <- function(df, Y, categoricalVars, 
                          trainSize=0.8, verbose=FALSE) {
  # parameters:
  #   Y - predictor variable
  #   categoricalVars - list of categorical variables to stratify on
  #   trainSize - percentage of df to put in training set
  #   verbose - boolean to togle printing the size of each set
  # returns:
  #   A list of a train data.frame and a test data.frame
  
  # add an id column (for error catching later)
  df$id <- as.numeric(rownames(df))
  
  # Train-Test Stratification
  sets = stratified(df, group=c(Y,categoricalVars), 
                    size=trainSize, bothSets=TRUE )
  
  # print out the size
  if(verbose == TRUE){
    print(paste("Train size:",nrow(sets[[1]])))
    print(paste("Test size:",nrow(sets[[2]])))
  }
  
  return(sets)
}
  
## --- fitEnsemble -----------------------------------------------------------------------

fitEnsemble <- function(train, test, Y, groupingVar, method='lasso', params=NULL) {
  # parameters:
  #   train - a data.frame to train the model
  #   test - a data.frame to test the model on 
  #   Y - outcome variable (column name) 
  #   groupingVar - name of the column used to identify groups
  #   method - 'lasso', 'tree', 'step'
  #   params - if method = 'tree' will be used as rpart.control
  # returns:
  #   list of info 
  
  l = length(levels(train[,groupingVar])) # number of groups
  groupNames = levels(train[,groupingVar]) # names of the groups
  
  AUC = list() # list to save testing AUC for each submodel
  train_AUC = list() # list to save training AUC for each submodel
  
  # create data structure to store predictions
  savedPreds = c()
  train_savedPreds = c()
  models = list()
  
  # iterating through the groups training/testing models
  for (j in 1:l){
    
    # progress bar
    if (j ==1) {print(paste(method, "ensemble"))}
    cat(paste0(" --> ", j, "/", l))
    if (j == l){ cat('\n') } # to end line
    
    # get training and testing data
    i <- groupNames[j]
    Xvars <- (colnames(train))[(colnames(train) != groupingVar)]
    trainSet <- train[which(train[,groupingVar]==i),Xvars]
    testSet <- test[which(test[,groupingVar]==i),Xvars]

    switch(method,
           'lasso' = {
             models[[i]] <- fitSimpleLasso(train=trainSet, test=testSet, 
                                           Y=Y, type.measure="auc", refit=TRUE)
             yhat <- as.vector(models[[i]]$refit.model$predY)
             train_yhat <- as.vector(models[[i]]$refit.model$train_predY)
             AUC[[i]] <- models[[i]]$refit.model$AUC
             train_AUC[[i]] <- models[[i]]$refit.model$train_AUC
           },
           'tree' = {
             if (!is.null(params)) { # if parameters are specified
               models[[i]] <- fitSimpleTree(trainSet, testSet, Y, control=params)
             } else {models[[i]] <- fitSimpleTree(trainSet, testSet, Y)}
             yhat <- as.vector(models[[i]]$predY)
             train_yhat <- as.vector(models[[i]]$train_predY)
             AUC[[i]] <- models[[i]]$AUC
             train_AUC[[i]] <- models[[i]]$train_AUC
           },
           'step' = {
             models[[i]] <- fitSimpleStep(train=trainSet, test=testSet, Y=Y, trace=FALSE)
             yhat <- as.vector(models[[i]]$predY)
             train_yhat <- as.vector(models[[i]]$train_predY)
             AUC[[i]] <- models[[i]]$AUC
             train_AUC[[i]] <- models[[i]]$train_AUC
           },
           stop('invalid method argument')
           )
    
    names(yhat) <- rownames(testSet) 
    names(train_yhat) <- rownames(trainSet) 
    savedPreds <- c(savedPreds,yhat)
    train_savedPreds <- c(train_savedPreds, train_yhat)
  }
  
  # order the predictions by their indices
  orderedPreds <- savedPreds[order(as.numeric(names(savedPreds)))]
  train_orderedPreds <- train_savedPreds[order(as.numeric(names(train_savedPreds)))] 
  
  # put everything together to return
  returninfo = list(AUC = AUC,
                    predY = orderedPreds,
                    trueY = test[,Y],
                    groupingVar = test[,groupingVar],
                    train_AUC = train_AUC,
                    train_predY = train_orderedPreds,
                    train_trueY = train[,Y],
                    train_groupingVar = train[,groupingVar],
                    models = models)
  
  return(returninfo)
}

## --- Model Fitting Functions -----------------------------------------------------------------------

fitSimpleModel <- function(train, test, Y){
  # parameters:
  #   train - data.frame
  #   test - data.frame
  #   Y - name of the outcome variable (assumed binary)
  # returns:
  #   list of info
  
  # train model
  mod = glm(as.formula(paste(Y,"~ .")),family="binomial", data=train)
  
  preds = predict(object = mod, newdata = test, type='response')
  
  returninfo = list(AUC = getROC(p=preds,y=test[,Y]),
                    predY = preds,
                    trueY = test[,Y],
                    train_AUC = getROC(predict(mod), train[,Y]),
                    train_predY = predict(mod),
                    train_trueY = train[,Y],
                    model = mod)
  
  return(returninfo)
}

fitSimpleTree <- function(train, test, Y, control=rpart.control(maxsurrogate=0,maxcompete = 0,cp=0.001)){
  # parameters:
  #   train - data.frame
  #   test - data.frame
  #   Y - name of the outcome variable (assumed binary)
  #   control - an rpart.control object
  # returns:
  #   list of info
  
  mod = rpart(as.formula(paste(Y,"~ .")), 
                data=train,
                control=control,
                method='class',
                parms=list(split='gini'))

  preds = predict(mod, test, na.action=na.omit)[,'1']
  returninfo = list(AUC = getROC(p=preds,y=test[,Y]),
                    predY = preds,
                    trueY = test[,Y],
                    train_AUC = getROC(p=predict(mod)[,'1'],y=train[,Y]),
                    train_predY = predict(mod)[,'1'],
                    train_trueY = train[,Y],
                    model = mod)
}
                           
fitSimpleStep <- function(train, test, Y, trace=FALSE){
  # parameters:
  #   train - data.frame
  #   test - data.frame
  #   Y - name of the outcome variable (assumed binary)
  #   trace - true/false passed to the trace option in step()
  # returns:
  #   list of info
  
  mod = glm(as.formula(paste(Y,"~ .")),family="binomial", data=train)
  s = step(mod, k=log(nrow(train)), trace = trace)
  
  preds = predict(object = s, newdata = test, type='response')
  
  returninfo = list(AUC = getROC(p=preds,y=test[,Y]),
                    predY = preds,
                    trueY = test[,Y],
                    train_AUC = getROC(p=predict(s),y=train[,Y]),
                    train_predY = predict(s),
                    train_trueY = train[,Y],
                    model = s
  )
  return(returninfo)
}

fitSimpleLasso <- function(train, test, Y, type.measure="auc", refit=TRUE, seed=4){
  # parameters:
  #   train - data.frame
  #   test - data.frame
  #   Y - name of the outcome variable (assumed binary)
  #   type.measure - passed to the type.measure option in cv.glmnet
  #   refit - true/false, whether to refit coefficients after lasso feature selection
  #   seed - integer
  # returns:
  #   list of info
  
  set.seed(seed) # influences how lasso split training data for CV
  
  # making model matrix
  trainX <- model.matrix(as.formula(paste(Y,"~ . + 0")), data=train)
  testX <- model.matrix(as.formula(paste(Y,"~ . + 0")), data=test)
  
  # training lasso regression
  fit.lasso <- glmnet(x=trainX, y=as.factor(train[,Y]), family="binomial")
  cv.lasso = cv.glmnet(x=trainX, y=as.factor(train[,Y]), family="binomial", 
                       type.measure = type.measure, nfolds=5)
  
  # making predictions
  lambda = cv.lasso$lambda.1se
  #if (length(getLassoVariables(coef(cv.lasso, s=lambda))) == 0){
  #  lambda = cv.lasso$lambda.min
  #}
  preds = as.vector(predict(fit.lasso, testX, s=lambda, type='response'))
  
  if (refit == TRUE){
    selected_vars = getLassoVariables(coef(cv.lasso, s=lambda))
    new_train = data.frame(cbind(train[,Y], trainX[,selected_vars]))
    new_test = data.frame(cbind(test[,Y], testX[,selected_vars]))
    colnames(new_train)[1] <- paste(Y)
    colnames(new_test)[1] <- paste(Y)
    refit.model = fitSimpleModel(new_train, new_test, Y=Y)
  } else { refit.model = NULL }
  
  returninfo = list(AUC = getROC(p=preds, y=test[,Y]),
                    trueY = test[,Y],
                    predY = preds,
                    lasso_coefs = coef(cv.lasso, s=lambda),
                    model = fit.lasso, 
                    cv.model = cv.lasso,
                    refit.model = refit.model
                    )
  
  return(returninfo)
}


getLassoVariables <- function(coef_matrix){
  # summary:
  #   helper function for fitSimpleLasso()
  # parameters:
  #   coef_matrix - a sparse matrix of variable names and estimated coefficients 
  #                 from a lasso model
  # returns: 
  #   a list of the variables estimated as non-zero by lasso for the chosen lambda
  lasso_vars = rownames(coef_matrix)[as.matrix(coef_matrix) != 0] 
  lasso_vars = lasso_vars[!(lasso_vars %in% "(Intercept)")]
  return(lasso_vars)
}

## --- CrossValidation -----------------------------------------------------------------------

repeatedSplits <- function(df, Y, groupingVar=NULL, splits, type, method, params=NULL){
  # parameters:
  #   df - data.frame object of data
  #   Y - the name of the Y variable e.g. 'readmitIn30'
  #   splits - a list of vectors, each vector has indices of test data for a split
  #   groupingVar - the name of a variable f
  # returns:
  #   a list of models and other results
  
  if (!(method %in% c('lasso','step','tree'))) {stop('must specify valid model method')}
  if (!(type %in% c('global','ensemble'))) {stop('must specify valid model type')}
  if (type=='ensemble' && is.null(groupingVar)) {stop('must specify groupingVar if ensemble')}
  
  # place to save output info
  models = list()
  AUC = list()
  train_AUC = list()
  
  # iterate over the folds
  for (i in 1:length(splits)){
    
    # progress bar
    cat(paste0(" --> ", 100*(i/length(splits)),"%"))
    if (i == length(split)){ cat('\n') } # to end line
    
    # get training and testing sets 
    training = df[-splits[[i]],]
    testing = df[splits[[i]],]
    
    switch(type,
           'global' = {
             
             switch(method,
                    'lasso' = {models[[i]] <- fitSimpleLasso(training, testing, Y)},
                    'step' = {models[[i]] <- fitSimpleStep(training, testing, Y)},
                    'tree' = {
                      if (!is.null(params)) { # if parameters are specified
                        models[[i]] <- fitSimpleTree(training, testing, Y, control=params)
                      } else {models[[i]] <- fitSimpleTree(training, testing, Y)}
                    }, 
                    stop('invalid method argument'))
             
             if (method == 'lasso') {
               AUC[[i]] <- models[[i]]$refit.model$AUC
               train_AUC[[i]] <- models[[i]]$refit.model$train_AUC
             } else {
               AUC[[i]] <- models[[i]]$AUC
               train_AUC[[i]] <- models[[i]]$train_AUC
             }
             
           },
           'ensemble' = {
             
             models[[i]] = fitEnsemble(train=training,test=testing, 
                                       Y=Y, groupingVar=groupingVar, method=method)
             AUC[[i]] <- models[[i]]$AUC
             train_AUC[[i]] <- models[[i]]$train_AUC
           },
           stop('invalid type argument')
    )
  }
  
  return(list(models=models,AUC=AUC, train_AUC=train_AUC, splits=split))
}


kCrossValidation <- function(df, Y, groupingVar=NULL, folds, type, method, params=NULL){
  # parameters:
  #   df - data.frame object of data
  #   Y - the name of the Y variable e.g. 'readmitIn30'
  #   folds - a list of numbers same length as df 
  #   groupingVar - the name of a variable f
  # returns:
  #   a list of models and other results
  
  if (!(method %in% c('lasso','step','tree'))) {stop('must specify valid model method')}
  if (!(type %in% c('global','ensemble'))) {stop('must specify valid model type')}
  if (type=='ensemble' && is.null(groupingVar)) {stop('must specify groupingVar if ensemble')}
  folds <- c(folds) # make sure the folds is a vector
  
  # place to save output info
  models = list()
  AUC = list()
  train_AUC = list()
  
  # iterate over the folds
  for (i in 1:length(unique(folds))){
    
    # progress bar
    cat(paste0(" --> ", 100*(i/length(unique(folds))),"%"))
    if (i == length(unique(folds))){ cat('\n') } # to end line

    # get training and testing sets 
    training = df[which(folds!=i),]
    testing = df[which(folds==i),]
    
    switch(type,
           'global' = {
             
             switch(method,
                    'lasso' = {models[[i]] <- fitSimpleLasso(training, testing, Y)},
                    'step' = {models[[i]] <- fitSimpleStep(training, testing, Y)},
                    'tree' = {
                      if (!is.null(params)) { # if parameters are specified
                        models[[i]] <- fitSimpleTree(training, testing, Y, control=params)
                      } else {models[[i]] <- fitSimpleTree(training, testing, Y)}
                    }, 
                    stop('invalid method argument'))
             
             if (method == 'lasso') {
               AUC[[i]] <- models[[i]]$refit.model$AUC
               train_AUC[[i]] <- models[[i]]$refit.model$train_AUC
             } else {
               AUC[[i]] <- models[[i]]$AUC
               train_AUC[[i]] <- models[[i]]$train_AUC
             }
             
           },
           'ensemble' = {
             
             models[[i]] = fitEnsemble(train=training,test=testing, 
                                      Y=Y, groupingVar=groupingVar, method=method)
             AUC[[i]] <- models[[i]]$AUC
             train_AUC[[i]] <- models[[i]]$train_AUC
           },
           stop('invalid type argument')
    )
  }
  
  return(list(models=models,AUC=AUC, train_AUC=train_AUC, folds=folds))
}

