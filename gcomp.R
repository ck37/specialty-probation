lib <- c("SL.gam","SL.polymars","SL.stepAIC","SL.glmnet.0","SL.glmnet.0.25","SL.glmnet.0.75","SL.glmnet.0.5","SL.glmnet.1","SL.earth")

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
 
# -0.1070508



txt <- W_arrest2
control <- W_arrest2
txt$treatment <- 1
control$treatment <- 0
newdata <- rbind(W_arrest2,txt,control)
fit2 <- SuperLearner(Y=data_arrest$any_arrest,X=W_arrest2,newX=newdata,family="binomial", SL.library=lib)
fit2
preds1 <- fit2$SL.predict[355:708]
preds0 <- fit2$SL.predict[709:1062]
mean(preds1-preds0)

# 







txt <- W_viol
control <- W_viol
txt$treatment <- 1
control$treatment <- 0
newdata <- rbind(txt,control)
fit.viol <- SuperLearner(Y=data_viol$any_viol,X=W_viol,newX=newdata,family="binomial", SL.library=lib)
fit.viol
preds1 <- fit.viol$SL.predict[1:291]
preds0 <- fit.viol$SL.predict[292:582]

mean(preds1-preds0)

# 0.0001




txt <- W_viol2
control <- W_viol2
txt$treatment <- 1
control$treatment <- 0
newdata <- rbind(txt,control)
fit.viol2 <- SuperLearner(Y=data_viol$any_viol,X=W_viol2,newX=newdata,family="binomial", SL.library=lib)
fit.viol2
preds1 <- fit.viol2$SL.predict[1:291]
preds0 <- fit.viol2$SL.predict[292:582]

mean(preds1-preds0)

#  0.001010674
