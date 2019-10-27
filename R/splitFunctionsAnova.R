
# 3 functions for implementing Anova using rpart's interface for user
# specified split functions. Consistent with rpart's C implementation of
# Anoval splitting.
#
# Adapted from rpart/tests/user_splits.R in rpart package

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
initAnova <- function(y, offset, parms, wt) {
  if (!is.null(offset)) y <- y-offset
  list(y=y, parms=parms, numresp=1, numy=1,
       summary= function(yval, dev, wt, ylevel, digits ) {
         paste("  mean=", format(signif(yval, digits)),
               ", MSE=" , format(signif(dev/wt, digits)),
               sep='')
       })
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
evalAnova <- function(y, wt, parms) {
  wmean <- sum(y*wt)/sum(wt)
  rss <- sum(wt*(y-wmean)^2)
  list(label= wmean, deviance=rss)
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
splitAnova <- function(y, wt, x, parms, continuous) {
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
    return( list(goodness=goodness, direction=sign(lmean)) )
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
    return(list(goodness = (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2),
                direction = ux[ord]))
  }
}

