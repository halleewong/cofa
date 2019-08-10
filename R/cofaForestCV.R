## --- Summary ----------------------------------------------------------------
# Hallee Wong
#
# Summary: functions for convieniently running cofrequency analysis
# only many forests of trees

# runs analyze forest k times and saves the frequency matrix, total matrix for
# each run
cofaForestCV <- function(k=10,
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




