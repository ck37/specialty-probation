
lib <- c("SL.DSA","SL.polymars","SL.stepAIC","SL.glmnet.0","SL.glmnet.0.25","SL.glmnet.0.75","SL.glmnet.0.5","SL.glmnet.1","SL.earth")

load("data/analysis-dataset.RData")

names(data)

data_viol <- data[!is.na(data$any_violence),]
data_arrest <- data[!is.na(data$any_arrest),]

X_viol = subset(data_viol, select=-c(studyid,treatment, any_violence, any_arrest))
X_viol <- X_viol[,1:31]
names(X_viol)

X_arrest = subset(data_arrest, select=-c(studyid, treatment,any_violence, any_arrest))
X_arrest <- X_arrest[,1:31]
names(X_arrest)

# Convert factors to column indicators.
W_viol = data.frame(model.matrix(~ . -1 , X_viol))

W_arrest = data.frame(model.matrix(~ . -1 , X_arrest))


# fit the model to predict P(A|W)

gfit <- SuperLearner(Y=data_arrest$treatment,X=W_arrest,family="binomial", SL.library=lib)

pred.g1W <- predict(gfit,type='response')[[1]]
pred.g0W <- 1-pred.g1W

gAW <- c()
length(gAW) <- nrow(data_arrest)
gAW[data_arrest$treatment==1] <- pred.g1W[data_arrest$treatment==1]
gAW[data_arrest$treatment==0] <- pred.g1W[data_arrest$treatment==0]

summary(gAW)

wt <- 1/gAW

# HT estimator

mean(as.numeric(data_arrest$treatment==1)*wt*data_arrest$any_arrest)/mean(as.numeric(data_arrest$treatment==1)*wt)-mean(as.numeric(data_arrest$treatment==0)*wt*data_arrest$any_arrest)/mean(as.numeric(data_arrest$treatment==0)*wt)

# -0.3439978



gfit_viol <- SuperLearner(Y=data_viol$treatment,X=W_viol,family="binomial", SL.library=lib)

pred.g1W <- predict(gfit_viol,type='response')[[1]]
pred.g0W <- 1-pred.g1W

gAW <- c()
length(gAW) <- nrow(data_viol)
gAW[data_viol$treatment==1] <- pred.g1W[data_viol$treatment==1]
gAW[data_viol$treatment==0] <- pred.g1W[data_viol$treatment==0]

summary(gAW)

wt <- 1/gAW

# HT estimator

mean(as.numeric(data_viol$treatment==1)*wt*data_viol$any_viol)/mean(as.numeric(data_viol$treatment==1)*wt)-mean(as.numeric(data_viol$treatment==0)*wt*data_viol$any_viol)/mean(as.numeric(data_viol$treatment==0)*wt)

# -0.1191004

