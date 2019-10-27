
# summary: retrieves relevant rows/columns from categorical split matrix
# parameters:
#   treefit - rpart.object
#   varName - string name of variable of concern
#   newNames - list of name of levels of <varName> e.g. levels(data$varName)
# returns:
#   matrix with each col corresp. to a level of <varName>
#   each row corresponding to a node split that used <varName>
getMatrix <- function(treefit, varName, levelNames){

  if (is.null(treefit$splits)){
    warning("cannot getMatrix of categorical splits from root tree")
    return(NULL)
  }
  if (length(unique(treefit$splits[,'count'])) != length(treefit$splits[,'count'])){
    warning("rpart.object$splits has splits with duplicate counts. May be printing surrogate/competing splits")
    # NOTE: if maxcompete=0 is not present in rpart.control then rpart.object$splits
    # will include alternative/surrogate splits!
  }

  # get integer matrix
  imat <- treefit$csplit
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


#' Counts the appearances of all categories in a single tree
#' @param mat a matrix returned by \code{getMatix()}
#' @return An integer vector containing the number of splits involving each category
getAllFreqofLevel  <- function(mat){
  vec <- colSums(!is.na(mat))
  names(vec) <- colnames(mat)
  return(vec)
}

#' Counts the number of times two categories split together versus split apart necessary to calculate the cofa statistic for a pair of levels
#'
#' This helper function is used by \code{makeCofaMat()} and \code{cofaForest()}
#'
#' @param i the index of a level
#' @param j the index of a level
#' @param mat a split matrix returned by \code{getMatrix()}
#' @return A list containing two objects:
#' \code{together} is the number of times level \code{i} and level \code{j} of the categorical variable of interest split in the same direction
#' \code{apart} is the number of times level \code{i} and level \code{j} split apart.

getCounts <- function(i,j, mat){
  #browser()
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
  return(list(together=together, apart=apart))

}


#' Calculates the cofa statistic for a pair of levels
#'
#' This helper function is used by \code{makeCofaMat()}
#'
#' @param i an integer
#' @param j an integer
#' @param mat a matrix of numeric values
#' @return Calculated cofa statistic for level \code{i} and level \code{j} of the categorical variable of interest using the tree split matrix \code{mat}. The cofa statistic is the probability that levels \code{i} and \level{j} split together at a node where the categorical variable of interest is used to split.
getCofa <- function(i,j, mat){

  # get together and apart counts
  lst <- getCounts(i,j,mat)

  # compute statistic
  if ((lst[['together']] + lst[['apart']]) == 0){
    return(NA)
  }else{
    return(lst[['together']]/(lst[['together']]+lst[['apart']]))
  }
}

#' Create a matrix of cofa statistics from the split matric of a single tree
#'
#' @param mat a matrix or vector.
#' @return A matrix where cell (i,j) contains the cofa statistic for level i and level j of the categorical variable of interest. The cofa statistic is the probability that levels i and j split together at a node where the categorical variable of interest is used to split.
makeCofaMat <- function(mat){

  if (class(mat) != 'matrix'){ # if mat is a vector
    d <- length(mat)
    names <- names(mat)
  } else{ # if the mat is a matrix
    d <- ncol(mat)
    names <- colnames(mat)
  }

  # make empty matrix
  freqMat <- matrix(nrow=d, ncol=d)

  lapply(1:d, function(i) {
    lapply(1:i, function(j) {
      prob <- getCofa(i,j,mat)
      freqMat[i,j] <- prob
      freqMat[j,i] <- prob
    })
  })

  colnames(freqMat) <- name
  rownames(freqMat) <- names

  return(freqMat)
}
