
require(rpart)

# Summary: Fits forest and returns frequency matrix
# params:
#   ntree - (integer) number of tree to grow
#   cvar - (string) name of categorical variable to study
#   data - (data.frame) data to subset
#   yvar - (string) name of y variable column
#   xvars - (list) names of x variables to use
#   method - (list) user defined functions or key word
#   control - (rpart.control) object to pass to control in rpart
#   parms - (list) list of parms to pass to rpart
#   seed - set.seed()
#   indices - list of lists to fit each tree to
# returns:
#   list(freqMat=<matrix>, trees=<list of rpart.obj>)
analyzeForest <- function(formula, 
                          variable, 
                          data,
                          method='class',
                          control=rpart.control(maxsurrogate=0,maxcompete = 0),
                          parms=list(split='gini'),
                        
){
  
  if ( ntree==0 & is.null(indices)) stop("Need to provide indices or ntree")
  
  # should check all indices are between 1:nrow(data)
  
  if ( !is.null(indices) ){
    # make sure ntree is consistent with length of indices
    ntree = length(indices)
  }
  
  # store info about categorical variable of interest
  levelNames = levels(data[,cvar])
  l = length(levelNames)
  
  # matrix to keep track of counts together, apart
  togetherMat <- matrix(0, nrow=l, ncol=l)
  apartMat <- matrix(0, nrow=l, ncol=l)
  colnames(togetherMat) = rownames(togetherMat) = levelNames
  colnames(apartMat) = rownames(apartMat) = levelNames
  
  treeList <- list()
  indexList <- list()
  
  for (k in 1:ntree){
    
    # get indices of data to use
    if (is.null(indices)){
      # re-sample data
      idx <- sample(1:nrow(data), size=nrow(data), replace=TRUE)
    } else {idx <- indices[[k]]}
    
    # fit tree
    tree.obj <- rpart(formula, data,
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
    localLevels <- levels(data[idx,cvar])
    
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
  
  return(list(freqMat = freqMat, # frequency matrix for individual objects
                trees = treeList, # list of tree objects
                totalMat = (togetherMat+apartMat), # matrix of count of times 2 levels used at node
                indices = indexList)) # list of lists of indices used to fit trees
  
}

