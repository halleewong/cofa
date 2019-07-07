# SplitFunctionsGiniRandom.R
# Hallee Wong
#
# Split function based on the gini split functions, considers only a subset of variables 
# at each node (sqrt(<number of variables>)) to split on (all other goodness is set to 0)
#

# The split function, called once per split variable per node
# 
# parameters:
#   y - vector of reponse values length n
#   wt - vector of weights
#   x - vector of x values for split variable being considered
#   params - vector of user parameters passed from rpart() call
#            **Should include k=<number of variables>**
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
#
splitGiniRandom <- function(y, wt, x, parms, continuous) {
  
  if (!('k' %in% names(parms))) stop ("parms = list(k=<number of variables>)")
  
  k = parms$k # number of variables 
  number <- sample(1:k,1)
  
  if (continuous) { 
    # continuous x variable
    n <- length(y)
    
    left.sum <- cumsum(y)[-n] # y1, y1+y2, y1+y2+y3,...
    
    # number example in each child
    left.n  <- cumsum(rep(1,n))[-n] # num examples going left
    right.n <- rep(n,n-1) - left.n # num examples going right
    
    # percent of examples of class 1
    left.p1 <- left.sum/left.n
    right.p1 <- (rep(sum(y),n-1) - left.sum)/right.n
    p1 <- sum(y) / length(y)  
    
    # gini (impurity metric)  
    I.A = rep((2 * p1 * (1 - p1)), n-1) # gini of parent node
    I.left = 2 * left.p1 * (1 - left.p1)
    I.right = 2 * right.p1 * (1 - right.p1)
    
    goodness <- ( length(y)*I.A - (left.n)*(I.left) - (right.n)*(I.right) )
    direction <- ifelse(left.p1 < right.p1, -1, +1)
    
    if (number > sqrt(k)){
      return(list(goodness=rep.int(0,length(goodness)), direction=direction))
    } else { 
      return(list(goodness = goodness, direction = direction))
    }
      
    
  }
  else {
    # categorical x variable
    ux <- sort(unique(x)) # list of group names
    nums <- tapply(rep(1,length(x)), x, sum) # num of ex. per group
    ysum  <- tapply(y, x, sum) # total y by group
    means <- ysum/nums # mean value per group, names are group #
    
    # For binary y, we can order the categories by their means
    ord <- order(means)
    n <- length(ord) # k the number of groups
    
    left.sum <- cumsum(ysum[ord])[-n] # sum for y in left child
    
    # number of examples per group
    left.n <- cumsum(nums[ord])[-n]
    right.n <- length(y) - left.n
    
    # percent of examples of class 1 
    # ... now reusing code from continuous example
    left.p1 <- left.sum/left.n
    right.p1 <- (rep(sum(y),n-1) - left.sum)/right.n
    p1 <- sum(y) / length(y) 
    
    # gini (impurity metric)  
    I.A = rep(( 2 * p1 * (1 - p1)), n-1) # gini of parent node
    I.left =  2 * left.p1 * (1 - left.p1)
    I.right =  2 * right.p1 * (1 - right.p1)
    
    goodness <- ( length(y)*I.A - (left.n)*(I.left) - (right.n)*(I.right) )
    
    if (number >= sqrt(k)){
      return(list(goodness=rep.int(0,length(goodness)), direction=ux[ord]))
    } else { 
    return(list(goodness=goodness, direction=ux[ord]))
    }
    
  }
}