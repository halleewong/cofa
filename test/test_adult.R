# Summary: test cofa functions using the adult dataset from UCI
# Date Created: 10/26/19
# Updates:

# Load packages ---------------------------------------------------------------

library(cofa)

library(data.table)

library(rpart.plot)
library(randomForest)

# Load Data --------------------------------------------------------------------

df <- fread("/data/adult.csv", stringsAsFactors = T)
colnames(df) <- c("age", "workclass", "fnlwgt", "education",
                  "education_num", "marital_status", "occupation",
                  "relationship", "race", "sex", "capital_gain",
                  "capital_loss", "hours_per_week",
                  "native_country", "income")

# Data cleaning ---------------------------------------------------------------

list0 <- c("age", "workclass", "education", "marital_status", "occupation",
           "relationship", "race", "sex", "hours_per_week", "native_country")

outcome <- 'income'

df0 <- df[ , c(outcome, list0), with=F]
df0[ , bin_income := ifelse(income==">50K", 1, 0)]

# Single CART tree -------------------------------------------------------------

thisFormula <- paste(outcome, '~', paste0(c(list0), collapse = " + "))

m <- rpart( thisFormula, data=data.frame(df0),
            control = rpart.control(cp=0.01, maxcompete = 0, maxsurrogate = 0),
            method = 'class')

rpart.plot(m)

# Single Random Forest ---------------------------------------------------------

rf <- randomForest::randomForest(f, data=data.frame(df0),importance=T)

importance(rf)

# Cofa ------------------------------------------------------------------------

col.category <- "education"
col.outcome <- "bin_income"

## Calculating cofa statistics on a single random forest

ptm <- proc.time()
test <- cofaForest(ntree=100,
                   cvar=col.category,
                   yvar=col.outcome,
                   xvars=list0,
                   data=df0)
proc.time() - ptm
# ~35sec for 10-tree random forest

vizCoFreqMat(test$freqMat, order=TRUE)

## Compare cofa stats to distribution of stats when labels are shuffled

rf_results <- cofaTest(k=100,
                       ntree=100,
                       cvar=col.category,
                       yvar=col.outcome,
                       xvars=list0,
                       data=df0)

n_categories <- length(unique(df0[[col.category]]))
n_tests <- (n_categories^2 - n_categories)/2
cutoff_w_bonf <- abs(qnorm(0.025/n_tests))

vizCoFreqMat(rf_results$primary$freqMat, order=TRUE)

vizMaskedMatrix(rf_results$primary,
                rf_results$nullhyp,
                metric="zScore", cutoff=cutoff_w_bonf, order=F)



