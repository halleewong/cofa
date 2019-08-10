## --- Summary ----------------------------------------------------------------
# Hallee Wong
#
# Summary: functions for convieniently running cofrequency analysis
# only many forests of trees
#
## --- Dependencies -----------------------------------------------------------
#source("./ccfAnalysis.R")

## --- New Functions ----------------------------------------------------------

# runs analyze forest k times and saves the frequency matrix, total matrix for 
# each run
runTrials <- function(k=10, 
                      dat=df1, 
                      xvars=c('c1','x1','x2'), 
                      control=rpart.control(cp=0.001,maxsurrogate=0,
                                            maxcompete=0,xval=0),
                      method='class',
                      parms=list(split='gini') ){
  
  results = list()
  totals = list()
  
  for (i in 1:k){
    print(paste("trial",i))
    temp <- analyzeForest(ntree=100, 
                          cvar='c1', yvar='y', xvars=xvars,
                          df=dat,
                          method=method, 
                          control=control,
                          parms=parms,
                          indices=NULL)
    results[[i]] <- temp$freqMat
    totals[[i]] <- temp$totalMat
  }
  
  return(list(results=results,totals=totals))
}

# takes the list of frequency matrices <results> and flattens them for 
# histogram plotting
gatherFMValues <- function(results,k=10){
  values = c()
  for (i in 1:length(results)){
    fm <- results[[i]]
    values <- c(values,as.numeric(fm[lower.tri(fm, diag=FALSE)]))
  }
  
  tbl <- data.frame(values=values)
  tbl$label <- as.factor(unlist(lapply(1:10, FUN=rep, length(values)/k ), 
                                recursive=TRUE))
  
  return(tbl)
} 

# plots a simple histogram of the aggregated fm values
plotFMValues <- function(df=tbl){
  ggplot(df, aes(x=values)) + xlim(-0.1,1.1) +
    geom_histogram(binwidth=0.01) +
    theme_minimal()
}

# plots density curves for each forest in the trial
plotTrials <- function(df=tbl){
  ggplot(df, aes(x=values, color=label)) + xlim(-0.1,1.1) +
    geom_density(na.rm=TRUE) +
    theme_minimal()
}


## --- New Functions ----------------------------------------------------------

# summary: Retrieves all test statistics for level1-level2 from the results
# parameters:
#   results - 
#   level1, level2 - 
# returns: vector of values
# 
valueDist <- function(results, level1, level2){
  vals = c()
  for (i in 1:length(results)){
    vals = c(vals, results[[i]]$fm[level1,level2])
  }
  return(vals)
}


# summary: makes historgram ggplot showing test statistics and distribution of 
#   values
# parameters:
#   result0 - list with fm and tot matrices of the results for th real data
#   results - list of list(fm=matrix, tot=matrix)
#   level1, level2 - names of levels (must match colnames in the matrices)
# returns: 
#   ggplot object
#
plotLevelsDist <- function(result0, results, level1, level2, binwidth=0.01, normal_distribution=TRUE){
  temp = data.frame(vals=valueDist(results, level1, level2))
  p <- ggplot(temp, aes(x=vals)) +
    geom_vline(xintercept=0.5, col='gray', lty=2) +
    geom_histogram(binwidth=binwidth) +
    geom_vline(xintercept=result0$fm[level1,level2], col=wes_palette('Zissou')[2], size=1) +
    scale_x_continuous(breaks=seq(0,1,0.2), limits=c(-0.01,1.01)) +
    scale_y_continuous(breaks=seq(0,70,20), limits=c(0,75)) +
    theme_light() + 
    labs(title=paste(level1, level2, 
                     ": stat = ", round(result0$fm[level1,level2],3),
                     ", total = ", result0$tot[level1,level2] ),
         x="Value", y="Count") +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          text = element_text(size=20))
    
  if (normal_distribution == TRUE){
      p <- p + stat_function(color=wes_palette('Zissou')[4], size=1, fun = function(x) dnorm(x, mean = mean(temp$vals), 
                                            sd = sd(temp$vals)) * nrow(temp) * binwidth)
  }

  return(p)
}

# summary: finds average statistic value for every pair of levels
# parameters
#   result - list of pairs of matrices
# returns: 
#   matrix of mean values
meanMatrix <- function(result){
  fmat = result[[1]]$fm
  for (i in 2:length(result)){
    fmat = fmat + result[[i]]$fm
  }
  return(fmat/length(result))
}

# summary: helper function for making various plots, where level names not important
# parameters:
#   result - list of pairs of matrices (fm and tot)
# returns:
#   data.frame object with all values from each null-hypothesis distribution 
aggregateDistributions <- function(result){
  
  levelNames = rownames(result[[1]]$fm) 
  allDist = c()
  group = c()
  level1name = c()
  level2name = c()
  
  count = 1
  for (m in 2:length(levelNames)){
    for (n in 1:(m-1) ){
      
      # pair of levels
      level1 = levelNames[n]
      level2 = levelNames[m]
      
      # get statistic value from all trials
      vals = c()
      for (i in 2:length(result)){
        vals = c(vals, result[[i]]$fm[level1,level2])
      }
      allDist = c(allDist, vals)
      group = c(group, rep(count, times=length(vals)))
      level1name = c(level1name, rep(level1, times=length(vals)))
      level2name = c(level2name, rep(level2, times=length(vals)))
      count = count + 1
    }
  }
  
  return(data.frame(values=allDist, label=factor(group), level1=level1name, level2=level2name))
}

## --- Scores -----------------------------------------------------------------

# summary: calculates a p-value by counting number of values in the 
# distribution that are more extreme than the value given
# parameters:
#   value - real number
#   dist - vector of values 
# returns: single value between 0 and 1
#
pValue <- function(value, dist){
  vals = dist - mean(dist)
  val0 = value - mean(dist)
  return((sum(vals <= -abs(val0)) + sum(vals >= abs(val0)))/length(dist))
}

# summary: calculates a z-score using the mean and sd of the dist
# parameters:
#   value - real number
#   dist - vector of values 
# returns: single value 
#
zScore <- function(value, dist){
  mean = mean(dist)
  sd = sd(dist)
  return((value - mean)/sd)
}


# summary: Calculates scores for all pairs of levels
# parameters:
#   result0 - a list with a fm matrix and tot matrix
#   results - a list of pairs of matrices
#   pValue, zScore - boolean for type of metric
# returns: a matrix of z-scores or p-values
#
metricMat <- function(result0, results, metric){

  if (!(metric %in% c("pValue","zScore"))){
    warning("Either pValue or zScore must be TRUE but not both")
    }

  # create empty matrix
  levels = colnames(result0$fm)
  mat = matrix(NA, nrow = length(levels), ncol=length(levels))
  colnames(mat) = rownames(mat) = levels
  
  # claculate scores
  for (i in 1:length(levels)){
    for (j in 1:i-1){
      level1 = levels[i]
      level2 = levels[j]
      
      if (metric == "zScore"){
        z = zScore(value=result0$fm[level1,level2],
                   dist=valueDist(results,level1,level2))
      } 
      else if (metric == "pValue"){
        z = pValue(value=result0$fm[level1,level2],
                   dist=valueDist(results,level1,level2))
      } else {z = NA}
      
      mat[level1,level2] = mat[level2,level1] = z
    }
  }
  
  return(mat)
}

# returns: 
# returns:
#   matrix object
getMaskedMat <- function(result0, results, metric, cutoff){
  
  if (metric=="pValue"){
    mask = abs(metricMat(result0, results, metric=metric)) < cutoff
  } else if (metric=="zScore"){
    mask = abs(metricMat(result0, results, metric=metric)) > cutoff
  }
  masked_mat = result0$fm
  masked_mat[!mask] <- NA
  diag(masked_mat) <- NA
  
  return(masked_mat)
}

# summary: Returns a 
# parameters:
#   result0 - list with fm and tot matrix
#   results - list of pairs of fm and tot matrices
#   metric - either "pValue" or "zScore"
#   cutoff - will be upper bound if using p-value and lower bound if using 
#            z-score
#   order - boolean passed to vizCoFreqMat
#   size - text size
# returns:
#   ggplot object of the matrix
vizMaskedMatrix <- function(result0, results, metric, cutoff, order, size=1){
  
  masked_mat = getMaskedMat(result0, results, metric, cutoff)
  
  if (order == FALSE){
    p <- vizCoFreqMat(mat, order=FALSE) + 
      geom_text(data = meltForViz(get_lower_tri(masked_mat)), 
                aes(label=signif(value,1)), 
                size=size, alpha=0.7) 
  } else {
    hc <- cluster_mat(result0$fm)
    masked_mat_ordered <- masked_mat[hc$order, hc$order]
    
    p <- vizCoFreqMat(result0$fm, order=TRUE) +  
      geom_text(data = meltForViz(get_lower_tri(masked_mat_ordered)), 
                aes(label=signif(value,1)), 
                size=size, alpha=0.7)  
  }
  return(p)
}