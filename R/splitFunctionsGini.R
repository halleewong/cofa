# SplitFunctionsGini.R
# Hallee Wong
#
# Set of 3 functions for implementing the gini split using rpart()'s interface for user
# specified split functions. Confirmed to produce the same splits at the C implementation
# of gini split for classification.
#
# TO DO
# - add print and text functions in the init function
#   to match rpart(method='class') output
#

# The init function
# parameters:
#   y - vector (?) of response values
#   params - specified by user in call of rpart(..., params=...)
#   wt - weight vector from the call 
# returns:
#   a list containing y, params, numresp, numy and some optional functions
#     y - vector of response values (possible updated if offset)
#     numresp - number of values produced by the eval routine's "label"    
#     numy - number of columns for y (should be 1)
#   optional functions
#     summary - (optional) function to produce 1-3 line summary for the node 
#               used by summary.rpart
#     print - (optional) function which will produce a one line summary used
#             by print
#     text - (optional) function 
#
initGini <- function(y, offset, parms, wt) {
  # offset parameter should be null
  if (is.null(offset)) offset <- 0
  if (any(y != 0 & y!= 1)) stop ('response must be 0/1')
  
  summaryFunction <- function (yval, dev, wt, ylevel, digits){
    
    #(paste("yval=", yval," dev=", dev, " wt=", wt, 
    #      " ylevel=", ylevel, " digits=", digits, sep=''))
    
    #yval is essentiall sum(y)/wt
    predGroup <- ifelse(yval < 0.5, 0, 1)
    
    nodeprob <- signif(wt / n, digits)
    
    # class counts
    class0 = paste(signif((1-yval) * wt, digits))
    class1 = paste(signif(yval * wt, digits))
    
    # class probabilities
    prob0 = paste(format(round(1-yval, 3), nsmall = 3))
    prob1 = paste(format(round(yval, 3), nsmall = 3))
    
    # padding for left justification
    l = max(length(class0),length(class1),5)
    
    fmt = paste0("%",l,"s")
    
    temp1 <- paste(sprintf(fmt, class0), sprintf(fmt, class1))
    temp2 <- paste(sprintf(fmt, prob0), sprintf(fmt, prob1))
    
    paste("  predicted class=", format(predGroup, justify = "left"), 
          "  expected loss=", signif(yval, digits), 
          "  P(node) =", format(nodeprob), 
          "\n    class counts: ", format(temp1, justify="right"), 
          "\n   probabilities: ", format(temp2, justify="right"), sep='') 
  }
  environment(summaryFunction) <- .GlobalEnv
  
  list(y=y, parms=parms, numresp=1, numy=1, 
       summary = summaryFunction)
}

# The 'evaluation' function, called once per node.
#
# parameters:
#   y - vector (?) of response values
#   params - specified by user in call of rpart(..., params=...)
#   wt - vector of weights input by user in the rpart() call 
# returns: 
#   a label (1 or more elements long) for labeling each node
#   and a deviance (of length 1) 
#     - equal to 0 if the node is "pure" in some sense (unsplittable)
#     - does not need to be a deviance: any measure that gets larger
#       as the node is less acceptable is fine.
#     - however, the measure underlies cost-complexity pruning
#
# for classification I think deviance should be # misclassified at node
#
evalGini <- function(y, wt, parms) {
  # probability y is 1 at node
  p1 <- sum(y)/length(y)
  
  if (p1 < 0.5){
    # node is classified as class 0
    d = sum(y == 1)
  } else {
    # node is classified as class 1
    d = sum(y == 0)
  }
  
  return(list(label = p1, deviance = d))
}

# The split function, called once per split variable per node
# 
# parameters:
#   y - vector of reponse values length n
#   wt - vector of weights
#   x - vector of x values for split variable being considered
#   params - vector of user parameters passed from rpart() call
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
splitGini <- function(y, wt, x, parms, continuous) {
  
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
    #print(goodness)
    
    return( list(goodness = goodness, 
                 direction = ifelse(left.p1 < right.p1, -1, +1)) )
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
    #print(goodness)
    
    return( list(goodness = c(goodness), direction = ux[ord])
    )
  }
}