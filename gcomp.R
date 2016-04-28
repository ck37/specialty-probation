
load("data/analysis-dataset.RData")

names(data)

data_viol <- data[!is.na(data$any_violence),]
data_arrest <- data[!is.na(data$any_arrest),]

X_viol = subset(data_viol, select=-c(treatment, studyid, any_violence, any_arrest))
X_viol <- X_viol[,1:31]
names(X_viol)

X_arrest = subset(data_arrest, select=-c(treatment, studyid, any_violence, any_arrest))
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

txt <- W
control <- W
txt$any_arrest==1
control$any_arrest==0

preds1 <- predict(fit,type='response',newdata=txt)
preds0 <- predict(fit,type='response',newdata=control)

mean(preds1-preds0)

