lib <- c("SL.DSA","SL.polymars","SL.stepAIC","SL.glmnet.0","SL.glmnet.0.25","SL.glmnet.0.75","SL.glmnet.0.5","SL.glmnet.1","SL.earth")

load("data/analysis-dataset.RData")

names(data)

data_viol <- data[!is.na(data$any_violence),]
data_arrest <- data[!is.na(data$any_arrest),]

X_viol = subset(data_viol, select=-c(studyid, any_violence, any_arrest))
X_viol <- X_viol[,1:31]
names(X_viol)

X_arrest = subset(data_arrest, select=-c(studyid, any_violence, any_arrest))
X_arrest <- X_arrest[,1:31]
names(X_arrest)

# Convert factors to column indicators.
W_viol = data.frame(model.matrix(~ . -1 , X_viol))

W_arrest = data.frame(model.matrix(~ . -1 , X_arrest))


# Can enable this to check for missing data, or run manually.
if (F) {
  # Review missing data in W.
  apply(W, MARGIN=2, FUN=function(col) { sum(is.na(col)) })
  
  # Review missing data in our treatment indicator.
  table(data$treatment, useNA="ifany")
  # Review missing data in our outcome indicator.
  table(data$any_arrest, useNA="ifany")
}



fit <- SuperLearner(Y=data_arrest$any_arrest,X=W_arrest,family="binomial", SL.library=lib)

txt <- W_arrest
control <- W_arrest
txt$treatment <- 1
control$treatment <- 0

preds1 <- predict(fit,type='response',newdata=txt)[[1]]
preds0 <- predict(fit,type='response',newdata=control)[[1]]

mean(preds1-preds0)

# -0.141365


fit.viol <- SuperLearner(Y=data_viol$any_viol,X=W_viol,family="binomial", SL.library=lib)

txt <- W_viol
control <- W_viol
txt$treatment <- 1
control$treatment <- 0

preds1 <- predict(fit.viol,type='response',newdata=txt)[[1]]
preds0 <- predict(fit.viol,type='response',newdata=control)[[1]]

mean(preds1-preds0)

# -0.006305265

