
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
  }
  
  return(mat)
}

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

