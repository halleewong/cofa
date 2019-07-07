## --- Summary ----------------------------------------------------------------
# Hallee Wong
#
# Summary: functions for implementing co-frequency analysis of the levels of
# categorical variables in trees, forests of trees and random forests of trees

## --- Dependencies -----------------------------------------------------------

# Required packages
require(ggplot2)
require(reshape2)
require(rpart)

# palette (for viz function)
require(wesanderson)
zpal = wes_palette("Zissou")

require(viridis)

## --- New Functions ----------------------------------------------------------

# summary: retrieves relevant rows/columns from categorical split matrix
# parameters:
#   treefit - rpart.object
#   varName - string name of variable of concern
#   newNames - list of name of levels of <varName> e.g. levels(data$varName)
# returns:
#   matrix with each col corresp. to a level of <varName>
#   each row corresponding to a node split that used <varName>
# NOTE: if maxcompete=0 is not present in rpart.control then
# rpart.object$splits will include alternative/surrogate splits!
#
getMatrix <- function(treefit, varName, levelNames){

  if (is.null(treefit$splits)){
    warning("cannot getMatrix of categorical splits from root tree")
    return(NULL)
  }
  if (length(unique(treefit$splits[,'count'])) != length(treefit$splits[,'count'])){
    warning("rpart.object$splits has splits with same counts. May be printing surrogate/competing splits")
    #print(treefit$splits)
  }

  # get integer matrix
  imat = treefit$csplit
  imat[imat == 2] <- NA

  # select rows and columns relevant to the specified variable
  nodes <- treefit$splits[which(rownames(treefit$splits) == varName),'index']

  if(length(nodes) < 1){
    print("no nodes with specified variable")
    return(NULL)
  }

  mat <- imat[nodes,1:length(levelNames)]

  if ( class(mat) == 'matrix' ){
    colnames(mat) <- levelNames
  } else {
    names(mat) <- levelNames
    #print(class(mat))
    #print(names(mat))
  }

  return(mat)
}

##----------------------------------------------------------------------------------------------
## Methods for individual trees
##----------------------------------------------------------------------------------------------

# summary: counts appearances of level in a single tree
# parameters:
#   i - (integer) index of level OR (string) name of level
#   mat - a matrix from getMatrix()
# returns:
#   (integer) number of splits involving the level
#   i.e. the number of rows with 1 or 3
getFreqofLevel <- function(mat,i){
  return( sum(!is.na(mat[,i])) )
}

# summary: counts appearances of all levels in a single tree
# parameters:
#   mat - a matrix from getMatrix()
# returns:
#   (integer vector) number of splits involving each level
#   i.e. number of rows with 1 or 3
getAllFreqofLevel  <- function(mat){
  vec <- colSums(!is.na(mat))
  names(vec) <- colnames(mat)
  return(vec)
}

# Summary: get the counts together and apart of two levels
#   from a matrix
# parmams:
#   i,j - (integer) indices or (string) names of levels to compare
#   mat - a labeled matrix or vector returned by getMatrix()
# returns:
#   list of the form list(togehter= ,apart= )
#   if the levels never appear together will return list(together=0, apart=0)
getCounts <- function(i,j, mat){


  if (class(mat) == 'matrix'){

    tbl = table(mat[,i] == mat[,j])

    if ("TRUE" %in% names(tbl)){
      together = tbl[["TRUE"]]
    } else { together = 0 }

    if ("FALSE" %in% names(tbl)){
      apart = tbl[["FALSE"]]
    } else { apart = 0 }

  } else {
    # if mat is a vector
    mat[is.na(mat)] <- 0 #otherwise will get an error

    if (mat[i] == mat[j]){
      together = 1
      apart = 0
    } else {
      together = 0
      apart = 1 }

  }
  return(list(together=together,apart=apart))

}


# Summary: Returns the probability of two levels splitting together if they
#   are at a node together (based on a single tree)
# parameters:
#   i,j - (integer) indices or (string) names of levels to compare
#   mat - a labeled matrix returned by getMatrix()
# returns:
#   the percent of splits that include the two levels in which they are
#   in the same group as opposed to different groups
getProbTogether <- function(i,j, mat){

  # get together and apart counts
  lst <- getCounts(i,j,mat)

  # compute statistic
  if ((lst[['together']] + lst[['apart']]) == 0){ return(NA) }
  else{
    return(lst[['together']]/(lst[['together']]+lst[['apart']]))
    }
}


# Summary: makes co-freq matrix from split matrix of single tree
# parameters:
#   mat - a matrix or vector returned by getMatrix()
# returns:
#   matrix of probability of together vs apart in a split
#   0 - if the levels appear in a split together they always branch apart
#   1 - if the levels appear together, they always branch together
#
makeCoFreqMat <- function(mat){

  if (class(mat) != 'matrix'){ # if mat is a vector
    d <- length(mat)
    names <- names(mat)
  } else{ # if the mat is a matrix
    d <- ncol(mat)
    names <- colnames(mat)
    }

  # make empty matrix
  freqMat <- matrix(nrow=d, ncol=d)

  for (i in 1:d){
    for (j in 1:i){
      freqMat[i,j] = freqMat[j,i] = getProbTogether(i,j,mat)
    }
  }

  colnames(freqMat) = rownames(freqMat) = names



  return(freqMat)
}

##----------------------------------------------------------------------------------------------
## Functions for forests of trees
##----------------------------------------------------------------------------------------------

# Summary: Fits forest and returns frequency matrix
# params:
#   ntree - (integer) number of tree to grow
#   cvar - (string) name of categorical variable to study
#   df - (data.frame) data to subset
#   yvar - (string) name of y variable column
#   xvars - (list) names of x variables to use
#   method - (list) user defined functions or key word
#   control - (rpart.control) object to pass to control in rpart
#   parms - (list) list of parms to pass to rpart
#   seed - set.seed()
#   indices - list of lists to fit each tree to
# returns:
#   list(freqMat=<matrix>, trees=<list of rpart.obj>)
analyzeForest <- function(ntree=100,
                          cvar, yvar, xvars,
                          df=cutData,
                          method='class',
                          control=rpart.control(maxsurrogate=0,maxcompete = 0),
                          parms=list(split='gini'),
                          indices=NULL,
                          verbose=FALSE
                          ){

  if ( ntree==0 & is.null(indices)) stop("Need to provide indices or ntree")

  # should check all indices are between 1:nrow(df)

  if ( !is.null(indices) ){
    # make sure ntree is consistent with length of indices
    ntree = length(indices)
    }

  # store info about categorical variable of interest
  levelNames = levels(df[,cvar])
  l = length(levelNames)

  # matrix to keep track of counts together, apart
  togetherMat <- matrix(0, nrow=l, ncol=l)
  apartMat <- matrix(0, nrow=l, ncol=l)
  colnames(togetherMat) = rownames(togetherMat) = levelNames
  colnames(apartMat) = rownames(apartMat) = levelNames

  treeList <- list()
  indexList <- list()

  for (k in 1:ntree){

    # progress bar
    if ( ntree > 10 & k %% round(ntree/10) == 0 ){
      cat( paste0(" --> ", 10*(k/round(ntree/10)),"%") )
      if (k == ntree){ cat('\n') } # to end line
      }

    # get indices of data to use
    if (is.null(indices)){
      # re-sample data
      idx <- sample(1:nrow(df), size=nrow(df), replace=TRUE)
    } else {idx <- indices[[k]]}

    # fit tree
    tree.obj <- rpart(formula=as.formula(paste(yvar,"~.")),
                      data=df[idx,c(yvar,xvars)],
                      control=control,
                      method=method,
                      parms=parms)

    # save information
    treeList[[k]] <- tree.obj
    indexList[[k]] <- idx

    if (is.null(tree.obj$splits)){
      print(paste("tree", k, "is root"))
      next
    }

    # in case a subset does not include all levels
    localLevels <- levels(df[idx,cvar])

    # get matrix of level splits
    splitMatrix <- getMatrix(tree.obj,varName=cvar,levelNames=localLevels)

    # if the cvar doesnt appear in any node
    if (is.null(splitMatrix)){
      print(paste("tree", k, "does not use cvar"))
      #print(tree.obj$splits)
      next }

    for (i in 1:length(localLevels)){
      for (j in 1:i){

        # get cofrequency counts
        lst <- getCounts(i,j,splitMatrix)
        together <- ifelse(is.na(lst[['together']]),0,lst[['together']])
        apart <- ifelse(is.na(lst[['apart']]),0,lst[['apart']])

        # update matrices
        name1 <- localLevels[[i]]
        name2 <- localLevels[[j]]
        togetherMat[name1,name2] = togetherMat[name2,name1] = togetherMat[name2,name1] + together
        apartMat[name1,name2] = apartMat[name2,name1] = apartMat[name2,name1] + apart
      }
    }

  }

  # calculate frequency matrix
  freqMat <- togetherMat / (togetherMat + apartMat)
  if (!isSymmetric(freqMat)) stop("Frequency matrix is not symmetric")

  if (verbose){
    # return list of info
    return(list(freqMat = freqMat, # frequency matrix for individual objects
                trees = treeList, # list of tree objects
                totalMat = (togetherMat+apartMat), # matrix of count of times 2 levels used at node
                indices = indexList)) # list of lists of indices used to fit trees
  } else {
    return(list(freqMat = freqMat,
                totalMat = (togetherMat+apartMat)))
  }

}

# Summary: prints all the trees from a analyzeForest result
# parameters:
#   result - (list) list include trees list of rpart.tree objects
# returns:
#   prints out all trees in forest
printAllTrees <- function(result){
  for (i in 1:length(result$trees)){
    rpart.plot(result$trees[[i]])
  }
}



