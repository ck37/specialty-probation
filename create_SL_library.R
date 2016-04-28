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


create.SL.glmnet <- function(alpha = c(0,0.25, 0.50, 0.75,1)) {
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet.', alpha[mm], '<- function(..., alpha = ', alpha[mm], ') SL.glmnet(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

create.SL.glmnet()

lib <- c("SL.DSA","SL.polymars","SL.stepAIC","SL.glmnet.0","SL.glmnet.0.25","SL.glmnet.0.75","SL.glmnet.0.5","SL.glmnet.1","SL.earth")
lib <- c("SL.polymars","SL.stepAIC","SL.glmnet.0","SL.glmnet.0.25","SL.glmnet.0.75","SL.glmnet.0.5","SL.glmnet.1","SL.earth")
