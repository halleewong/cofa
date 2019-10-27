
#' Fits a Random Forest and calculated CoFA statistics for all pairs of categories.
#'
#' @param ntree integer number of trees to fit within the random forest
#' @param yvar name of the column containing a binary outcome variable
#' @param cvar name of the
#' @param xvars vector of names of columns to use as predictors
#' @param data a data.frame including columns \code{yvar}, \code{cvar} and \code{xvars}
#' @param indices a list of vectors containing the indices of rows to be used when fitting each tree. If \code{NULL}, rows will be randomly sampled with replacement.
#'
#' @return A list of outputs including a matrix of CoFA statistics (\code{freqMat}), a matrix with the total number of times each pair of categories are used for splitting (\code{totalMat}), a list of rpart tree objects (\code{tree}) and a list of the indices used to fit each tree (\code{indices}).
#

cofaForest <- function(ntree=100,
                       cvar,
                       yvar,
                       xvars,
                       data,
                       control.cp = 0.001,
                       seed = 4,
                       indices=NULL
){

  df <- data.frame(data)

  lapply(colnames(df), function(col) {
    if ( ! col %in% colnames(df) ) stop(paste(col, "is not a column in data"))
  })

  if ( ntree==0 & is.null(indices)) stop("Need to provide indices or ntree")
  if ( is.null(indices) & !is.null(seed) ) set.seed(seed)

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
                      data=df[idx,unique(c(yvar,xvars,cvar))],
                      control=rpart.control(cp=control.cp, maxsurrogate=0, maxcompete = 0),
                      method=list(eval=evalGini, split=splitGiniRandom, init=initGini),
                      parms=list(k=length(unique(c(xvars,cvar)))))

    # save information
    treeList[[k]] <- tree.obj
    indexList[[k]] <- idx

    if (is.null(tree.obj$splits)){
      print(paste("Tree", k, "is a root"))
      next
    }

    # in case a subset does not include all levels
    localLevels <- levels(df[idx,cvar])

    # get matrix of level splits
    splitMatrix <- getMatrix(tree.obj,varName=cvar,levelNames=localLevels)

    # if the cvar doesnt appear in any node
    if (is.null(splitMatrix)){
      print(paste("Tree", k, "does not use", cvar))
      #print(tree.obj$splits)
      next
      }

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

  # return list of info
  return(list(freqMat = freqMat, # frequency matrix for individual objects
              trees = treeList, # list of tree objects
              totalMat = (togetherMat+apartMat), # matrix of count of times 2 levels used at node
              indices = indexList)) # list of lists of indices used to fit trees

}

