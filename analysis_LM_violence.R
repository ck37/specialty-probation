#VIOLENCE!!!!
library("SuperLearner")
library("tmle")

load("data/analysis-dataset.RData")
load("data/covariates.RData")

SL.library = c("SL.glm", "SL.step", "SL.glm.interaction", "SL.step.interaction")

#Dealing with missing data for outcome
data = data[!is.na(data$any_violence),]

#Data frame X with covariates and exposure
X = data.frame(A = data$treatment, data[,-c(grep("miss", names(data)), which(names(data)=="any_arrest"))])
#Data frames with A = 1 and A = 0
X1 = X0 = X
X1$A = 1
X0$A = 0
newdata = rbind(X, X1, X0)

#Superlearn for Qbar
Qinit = SuperLearner(Y = data$any_violence, X = X, newX = newdata, SL.library = SL.library, family = "binomial", )
n = dim(data)[1]
QbarAW = Qinit$SL.predict[1:n]
Qbar1W = Qinit$SL.predict[(n+1):(2*n)]
Qbar0W = Qinit$SL.predict[(2*n+1):(3*n)]

###SUBSTITUTION ESTIMATOR###
PsiHat.SS = mean(Qbar1W - Qbar0W)

#Data with only covariates and variables
W = subset(X, select = -c(A))

#Superlearn for g
gHatSL = SuperLearner(Y = data$treatment, X = W, SL.library = SL.library, family = "binomial")

#Pred prob of being experienced, given covariates
gHat1W = gHatSL$SL.predict
gHat0W = 1 - gHat1W

#Pred prob of obs experience, given covariates
gHatAW = rep(NA, n)
gHatAW[data$treatment == 1] = gHat1W[data$treatment == 1]
gHatAW[data$treatment == 0] = gHat0W[data$treatment == 0]

###IPTW ESTIMATOR###
PsiHat.IPTW = mean(as.numeric(data$treatment==1)*data$any_violence/gHatAW) - mean(as.numeric(data$treatment==0)*data$any_violence/gHatAW)

#Clever covariate for each subject
H.AW = (2*data$treatment-1)/ gHatAW

#Clever covariates at A = 1 and A = 0 for all subjects
H.1W = 1/gHat1W
H.0W = -1/gHat0W

#Update initial estiamte of Qbar using logistic regression Y ~ H.AW
logitUpdate = glm(data$any_violence ~ -1 +offset(qlogis(QbarAW)) + H.AW, family='binomial') 
eps = logitUpdate$coef

#Calculate predicted values for each subject under each treatment
QbarAW.star = plogis(qlogis(QbarAW)+ eps*H.AW) 
Qbar1W.star = plogis(qlogis(Qbar1W)+ eps*H.1W) 
Qbar0W.star = plogis(qlogis(Qbar0W)+ eps*H.0W)

###TMLE ESTIMATE###
PsiHat.TMLE = mean(Qbar1W.star) - mean(Qbar0W.star)

#Using function
PsiHat.TMLE2 = tmle(Y=data$any_violence, A=data$treatment, W=W, Q=cbind(Qbar0W, Qbar1W), g1W=gHat1W, family="binomial")
