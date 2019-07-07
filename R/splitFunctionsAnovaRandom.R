# SplitFunctionsAnovaRandom.R
# Hallee Wong
# 
# Split function based on the anova split functions, considers only a subset of variables 
# at each node (sqrt(<number of variables>)) to split on (all other goodness is set to 0)
#

# The split function, called once per split variable per node
# 
# parameters:
#   y - vector of reponse values length n
#   wt - vector of weights
#   x - vector of x values for split variable being considered
#   params - vector of user parameters passed from rpart() call
#           **Should include k=<number of variables>**
#   continuous - true/false 
# returns:
#   If continuous...
#   - the x vector is ordered, y vector is sorted to correspond (no missing)
#   - should return 2 vectors length (n-1) 
#      goodness - goodness of the split, larger numbers are better.
#                 0 means couldn't find any worthwhile split
#                 the ith value of goodness evaluates splitting 
#                 observations 1:i vs (i+1):n
#      direction - (-1) = send "y < cutpoint" to the left side of the tree
#                   (1) = send "y < cutpoint" to the right
# 
#   If categorical...
#   - x is a set of integers defining the groups for an unordered predictor
#   - should return a vector length k (# groups) and (k-1)
#       direction - gives the order to line the groups up in (by y mean)
#           so that only m-1 splits need to be evaluated rather than 2^(m-1)
#       goodness - vector of m-1 values
#
# Note: this is not a big deal, but making larger "mean y's" move towards
#   the right of the tree, as we do here, seems to make it easier to read.
#   
#   The reason for returning a vector of goodness is that the C routine
#   enforces the "minbucket" constraint. It selects the best return value
#   that is not too close to an edge.

splitAnovaRandom <- function(y, wt, x, parms, continuous) {
  
  if (!('k' %in% names(parms))) stop ("parms = list(k=<number of variables>)")

  k = parms$k # number of variables 
  number <- sample(1:k,1)
  
  # Center y
  n <- length(y)
  y <- y- sum(y*wt)/sum(wt)
  
  if (continuous) {
    # continuous x variable
    temp <- cumsum(y*wt)[-n]
    
    left.wt  <- cumsum(wt)[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp/left.wt
    rmean <- -temp/right.wt
    goodness <- (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2)
    
    if (number >= sqrt(k)){
      return(list(goodness=rep.int(0,length(goodness)), direction=sign(lmean)))
    } else {
      return( list(goodness=goodness, direction=sign(lmean)) )
    }
    
  }
  else {
    # Categorical X variable
    ux <- sort(unique(x))
    wtsum <- tapply(wt, x, sum)
    ysum  <- tapply(y*wt, x, sum)
    means <- ysum/wtsum
    
    # For anova splits, we can order the categories by their means
    #  then use the same code as for a non-categorical
    ord <- order(means)
    n <- length(ord)
    temp <- cumsum(ysum[ord])[-n]
    left.wt  <- cumsum(wtsum[ord])[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp/left.wt
    rmean <- -temp/right.wt
    goodness = (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2)

    if (number > sqrt(10)){
      return(list(goodness=rep.int(0,length(goodness)), direction=ux[ord]))
    } else {
      return(list(goodness = goodness, direction = ux[ord]))
    }
  }
  
}
