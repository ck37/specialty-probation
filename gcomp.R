
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

fit <- SuperLearner(Y=data$any_arrest,X=W,family="binomial", SL.library=lib)

txt <- W
control <- W
txt$any_arrest==1
control$any_arrest==0

preds1 <- predict(fit,type='response',newdata=txt)
preds0 <- predict(fit,type='response',newdata=control)

mean(preds1-preds0)

