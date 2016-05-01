lib <- c("SL.polymars","SL.stepAIC","SL.glmnet.0","SL.glmnet.0.25","SL.glmnet.0.75","SL.glmnet.0.5","SL.glmnet.1","SL.earth")

load("data/analysis-dataset.RData")

names(data)

data_viol <- data[!is.na(data$any_violence),]
data_arrest <- data[!is.na(data$any_arrest),]

X_viol = subset(data_viol, select=-c(studyid, any_violence, any_arrest))
X_viol <- X_viol[,1:30]

X_viol2 = subset(data_viol, select=-c(studyid, any_violence, any_arrest,PAI_somatization,PAI_schizophrenia,PAI_paranoia,PAI_depression))
X_viol2 <- X_viol2[,1:26]


X_arrest = subset(data_arrest, select=-c(studyid, any_violence, any_arrest))
X_arrest <- X_arrest[,1:30]

X_arrest2 = subset(data_arrest, select=-c(studyid,any_violence, any_arrest,PAI_somatization,PAI_schizophrenia,PAI_paranoia,PAI_depression))
X_arrest2 <- X_arrest2[,1:26]


# Convert factors to column indicators.
W_viol = data.frame(model.matrix(~ . -1 , X_viol))
W_viol2 = data.frame(model.matrix(~ . -1 , X_viol2))

W_arrest = data.frame(model.matrix(~ . -1 , X_arrest))
W_arrest2 = data.frame(model.matrix(~ . -1 , X_arrest2))


library(SuperLearner)

fit <- SuperLearner(Y=data_arrest$any_arrest,X=W_arrest,family="binomial", SL.library=lib)
fit
txt <- W_arrest
control <- W_arrest
txt$treatment <- 1
control$treatment <- 0


preds1 <- predict(fit,type='response',newdata=txt)[[1]]
preds0 <- predict(fit,type='response',newdata=control)[[1]]

mean(preds1-preds0)

# -0.04334923




fit2 <- SuperLearner(Y=data_arrest$any_arrest,X=W_arrest2,family="binomial", SL.library=lib)
fit2
txt <- W_arrest2
control <- W_arrest2
txt$treatment <- 1
control$treatment <- 0


preds1 <- predict(fit2,type='response',newdata=txt)[[1]]
preds0 <- predict(fit2,type='response',newdata=control)[[1]]

mean(preds1-preds0)

# -0.07241793








fit.viol <- SuperLearner(Y=data_viol$any_viol,X=W_viol,family="binomial", SL.library=lib)
fit.viol
txt <- W_viol
control <- W_viol
txt$treatment <- 1
control$treatment <- 0

preds1 <- predict(fit.viol,type='response',newdata=txt)[[1]]
preds0 <- predict(fit.viol,type='response',newdata=control)[[1]]

mean(preds1-preds0)

# -0.003660523





fit.viol2 <- SuperLearner(Y=data_viol$any_viol,X=W_viol2,family="binomial", SL.library=lib)
fit.viol2
txt <- W_viol2
control <- W_viol2
txt$treatment <- 1
control$treatment <- 0

preds1 <- predict(fit.viol2,type='response',newdata=txt)[[1]]
preds0 <- predict(fit.viol2,type='response',newdata=control)[[1]]

mean(preds1-preds0)

#  -0.00461374
