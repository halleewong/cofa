
#' Calculates CoFA statistics for a random forest on the provided data and simulated the null hypothesis distribution by calculating CoFA statistics on data with the categorical variable shuffled
#' @param k number of iterations to simulate null hypothesis CoFA statistics
#' @param data data.frame containing columns \code{xvars}, \code{yvar}, \code{cvar}
#' @param yvar name of column containing the outcome variable
#' @param cvar name of column containing the categorical variable of interest
#' @param xvar vector of other columns to use in fitting random forests
#' @param ntree number of trees per random forest
#' @param seed a random seed
#'
cofaTest <- function(k=10,
                     data,
                     xvars,
                     yvar,
                     cvar,
                     ntree=100,
                     seed=NULL){

  df <- data.frame(data)

  if ( !is.null(seed) ) set.seed(seed)

  original_labels <- df[[col.category]]

  # lists for saving the labels
  labels <- append(
    # first entry is the original label order
    list(original_labels),
    # reshuffle labels k times
    lapply(1:k, function(x){sample(original_labels, size = length(original_labels), replace=FALSE)} )
  )

  primary_forest <- cofaForest(ntree=ntree,
                               cvar=cvar,
                               yvar=yvar,
                               xvars=xvars,
                               data=data,
                               indices=NULL,
                               seed=NULL)

  fit_forest <- function(i){

    print(paste("Forest",i))

    # Temporarily change category label column
    df[[cvar]] <- labels[[i]]

    # fit a forest
    temp <- cofaForest(ntree=ntree,
                       cvar=cvar,
                       yvar=yvar,
                       xvars=xvars,
                       data=data,
                       indices=NULL,
                       seed=NULL)

    # save matrices to lists
    list(freqMat=temp$freqMat, totalMat=temp$totalMat)
  }

  nullhyp_results <- lapply(2:(k+1), FUN = fit_forest)
  #results <- pbmclapply(X=2:(k+1), FUN=fit_forest, mc.cores=6)

  return(list(primary=primary_forest,
              nullhyp=nullhyp_results,
              labels=labels))
}




