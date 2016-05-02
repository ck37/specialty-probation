
lib <- c("SL.polymars","SL.stepAIC","SL.glmnet.0","SL.glmnet.0.25","SL.glmnet.0.75","SL.glmnet.0.5","SL.glmnet.1","SL.earth")

load("data/analysis-dataset.RData")

names(data)

data_viol <- data[!is.na(data$any_violence),]
data_arrest <- data[!is.na(data$any_arrest),]

X_viol = subset(data_viol, select=-c(studyid,treatment, any_violence, any_arrest))
X_viol <- X_viol[,1:29]
names(X_viol)
X_viol2 = subset(data_viol, select=-c(studyid,treatment, any_violence, any_arrest,PAI_somatization,PAI_schizophrenia,PAI_paranoia,PAI_depression))
X_viol2 <- X_viol2[,1:25]


X_arrest = subset(data_arrest, select=-c(studyid, treatment,any_violence, any_arrest))
X_arrest <- X_arrest[,1:29]
names(X_arrest)
X_arrest2 = subset(data_arrest, select=-c(studyid, treatment,any_violence, any_arrest,PAI_somatization,PAI_schizophrenia,PAI_paranoia,PAI_depression))
X_arrest2 <- X_arrest2[,1:25]


# Convert factors to column indicators.
W_viol = data.frame(model.matrix(~ . -1 , X_viol))
W_viol2 <- data.frame(model.matrix(~ . -1 , X_viol2))

W_arrest = data.frame(model.matrix(~ . -1 , X_arrest))
W_arrest2 = data.frame(model.matrix(~ . -1 , X_arrest2))


# fit the model to predict P(A|W)

gfit <- SuperLearner(Y=data_arrest$treatment,X=W_arrest,family="binomial", SL.library=lib)
pred.g1W <- predict(gfit,type='response')[[1]]
pred.g0W <- 1-pred.g1W
gAW <- c()
length(gAW) <- nrow(data_arrest)
gAW[data_arrest$treatment==1] <- pred.g1W[data_arrest$treatment==1]
gAW[data_arrest$treatment==0] <- pred.g0W[data_arrest$treatment==0]
summary(gAW)
wt <- 1/gAW
mean(as.numeric(data_arrest$treatment==1)*wt*data_arrest$any_arrest)/mean(as.numeric(data_arrest$treatment==1)*wt)-mean(as.numeric(data_arrest$treatment==0)*wt*data_arrest$any_arrest)/mean(as.numeric(data_arrest$treatment==0)*wt)
# -0.2079611

mean(as.numeric(data_arrest$treatment==1)*wt*data_arrest$any_arrest)-mean(as.numeric(data_arrest$treatment==0)*wt*data_arrest$any_arrest)
# 



gfit2 <- SuperLearner(Y=data_arrest$treatment,X=W_arrest2,family="binomial", SL.library=lib)
pred.g1W <- predict(gfit2,type='response')[[1]]
pred.g0W <- 1-pred.g1W
gAW <- c()
length(gAW) <- nrow(data_arrest)
gAW[data_arrest$treatment==1] <- pred.g1W[data_arrest$treatment==1]
gAW[data_arrest$treatment==0] <- pred.g0W[data_arrest$treatment==0]
summary(gAW)
wt <- 1/gAW
mean(as.numeric(data_arrest$treatment==1)*wt*data_arrest$any_arrest)/mean(as.numeric(data_arrest$treatment==1)*wt)-mean(as.numeric(data_arrest$treatment==0)*wt*data_arrest$any_arrest)/mean(as.numeric(data_arrest$treatment==0)*wt)
# -0.2069941

mean(as.numeric(data_arrest$treatment==1)*wt*data_arrest$any_arrest)-mean(as.numeric(data_arrest$treatment==0)*wt*data_arrest$any_arrest)









gfit_viol <- SuperLearner(Y=data_viol$treatment,X=W_viol,family="binomial", SL.library=lib)
pred.g1W <- predict(gfit_viol,type='response')[[1]]
pred.g0W <- 1-pred.g1W
gAW <- c()
length(gAW) <- nrow(data_viol)
gAW[data_viol$treatment==1] <- pred.g1W[data_viol$treatment==1]
gAW[data_viol$treatment==0] <- pred.g0W[data_viol$treatment==0]
summary(gAW)
wt <- 1/gAW
mean(as.numeric(data_viol$treatment==1)*wt*data_viol$any_viol)/mean(as.numeric(data_viol$treatment==1)*wt)-mean(as.numeric(data_viol$treatment==0)*wt*data_viol$any_viol)/mean(as.numeric(data_viol$treatment==0)*wt)
# -0.009732972

mean(as.numeric(data_viol$treatment==1)*wt*data_viol$any_viol)-mean(as.numeric(data_viol$treatment==0)*wt*data_viol$any_viol)
# 0.01462514




gfit_viol2 <- SuperLearner(Y=data_viol$treatment,X=W_viol2,family="binomial", SL.library=lib)
pred.g1W <- predict(gfit_viol2,type='response')[[1]]
pred.g0W <- 1-pred.g1W
gAW <- c()
length(gAW) <- nrow(data_viol)
gAW[data_viol$treatment==1] <- pred.g1W[data_viol$treatment==1]
gAW[data_viol$treatment==0] <- pred.g0W[data_viol$treatment==0]
summary(gAW)
wt <- 1/gAW
mean(as.numeric(data_viol$treatment==1)*wt*data_viol$any_viol)/mean(as.numeric(data_viol$treatment==1)*wt)-mean(as.numeric(data_viol$treatment==0)*wt*data_viol$any_viol)/mean(as.numeric(data_viol$treatment==0)*wt)
# -0.01241035
mean(as.numeric(data_viol$treatment==1)*wt*data_viol$any_viol)-mean(as.numeric(data_viol$treatment==0)*wt*data_viol$any_viol)
# 0.0140722


