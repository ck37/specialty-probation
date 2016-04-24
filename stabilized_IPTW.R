
load("data/analysis-dataset.RData")

names(data)

# Handle missing outcome data.
# We might (arguably) assume that missing data means that no bad thing happened.
# There is probably some fancy missing data approach to use.
data$any_violence[is.na(data$any_violence)] = 0
data$any_arrest[is.na(data$any_arrest)] = 0

# Remove outcomes, assignment, and study id from dataset when creating X dataframe.
X = subset(data, select=-c(treatment, studyid, any_violence))
names(X)

# Convert factors to column indicators.
W = model.matrix(~ . -1 , X)
colnames(W)

# Can enable this to check for missing data, or run manually.
if (F) {
  # Review missing data in W.
  apply(W, MARGIN=2, FUN=function(col) { sum(is.na(col)) })
  
  # Review missing data in our treatment indicator.
  table(data$treatment, useNA="ifany")
  # Review missing data in our outcome indicator.
  table(data$any_arrest, useNA="ifany")
}


### SUPERLEARNER LIBRARY

library(SuperLearner)

SL.DSA <- function(Y, X, newX, family, obsWeights, maxsize = 2*ncol(X), maxorderint = 2, maxsumofpow = 2, Dmove = TRUE, Smove = TRUE, vfold = 5, ...) {
  require('DSA')
  dsaweights <- matrix(obsWeights, nrow = (vfold +1), ncol = nrow(X), byrow = TRUE)
  fit.DSA <- DSA(Y ~ 1, data = data.frame(Y, X), family = family, maxsize = maxsize, maxorderint = maxorderint, maxsumofpow = maxsumofpow, Dmove = Dmove, Smove = Smove, vfold = vfold, weights = dsaweights)
  pred <- predict(fit.DSA, newdata = newX)
  if(family$family == "binomial") { pred <- 1 / (1 + exp(-pred))}
  fit <- list(object = fit.DSA)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.DSA")
  return(out)
}

# 
predict.SL.DSA <- function(object, newdata, family, ...) {
  require('DSA')
  pred <- predict(object = object$object, newdata = newdata)
  if(family$family == "binomial"){ pred <- 1 / (1 + exp(-pred))}
  return(pred)
}

lib <- c("SL.DSA","SL.polymars","SL.glm","SL.randomForest","SL.step.forward","SL.step.interaction","SL.glmnet","SL.knn")

# fit the model to predict P(A|W)

gfit <- SuperLearner(Y=data$treatment,X=W,family="binomial", SL.library=lib)

pred.g1W <- predict(gfit,type='response')
pred.g0W <- predict(gfit,type='response')

gAW <- c()
length(gAW) <- nrow(data)
gAW[data$treatment==1] <- pred.g1W[data$treatment==1]
gAW[data$treatment==0] <- pred.g1W[data$treatment==0]

summary(gAW)

wt <- 1/gAW

# HT estimator

mean(as.numeric(data$treatment==1)*wt*data$any_arrest)/mean(as.numeric(data$treatment==1)*wt) 
-mean(as.numeric(data$treatment==0)*wt*data$any_arrest)/mean(as.numeric(data$treatment==0)*wt)




