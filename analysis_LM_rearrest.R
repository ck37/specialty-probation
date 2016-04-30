#ARREST!!!!
library("SuperLearner")
library("tmle")

load("data/analysis-dataset.RData")
load("data/covariates.RData")

SL.library = c("SL.step")

#Dealing with missing data for outcome
data = data[!is.na(data$any_arrest),]

#Data frame X with covariates and exposure
X = subset(data, select = -c(studyid, any_arrest, any_violence))
X = data.frame(model.matrix(~ . -1 , X))

#Data frames with A = 1 and A = 0
X1 = X0 = X
X1$treatment = 1
X0$treatment = 0
newdata = rbind(X, X1, X0)

#Superlearn for Qbar
Qinit = SuperLearner(Y = data$any_arrest, X = X, newX = newdata, SL.library = SL.library, family = "binomial", )
n = dim(data)[1]
QbarAW = Qinit$SL.predict[1:n]
Qbar1W = Qinit$SL.predict[(n+1):(2*n)]
Qbar0W = Qinit$SL.predict[(2*n+1):(3*n)]

###SUBSTITUTION ESTIMATOR###
PsiHat.SS = mean(Qbar1W - Qbar0W)

#Data with only covariates and variables
W = subset(X, select = -c(treatment))

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
PsiHat.IPTW = mean(as.numeric(data$treatment==1)*data$any_arrest/gHatAW) - mean(as.numeric(data$treatment==0)*data$any_arrest/gHatAW)

#Clever covariate for each subject
H.AW = (2*data$treatment-1)/ gHatAW

#Clever covariates at A = 1 and A = 0 for all subjects
H.1W = 1/gHat1W
H.0W = -1/gHat0W

#Update initial estiamte of Qbar using logistic regression Y ~ H.AW
logitUpdate = glm(data$any_arrest ~ -1 +offset(qlogis(QbarAW)) + H.AW, family='binomial') 
eps = logitUpdate$coef

#Calculate predicted values for each subject under each treatment
QbarAW.star = plogis(qlogis(QbarAW)+ eps*H.AW) 
Qbar1W.star = plogis(qlogis(Qbar1W)+ eps*H.1W) 
Qbar0W.star = plogis(qlogis(Qbar0W)+ eps*H.0W)

###TMLE ESTIMATE###
PsiHat.TMLE = mean(Qbar1W.star) - mean(Qbar0W.star)

#Using function
PsiHat.TMLE2 = tmle(Y=data$any_arrest, A=data$treatment, W=W, Q=cbind(Qbar0W, Qbar1W), g1W=gHat1W, family="binomial")
PsiHat.TMLE2$estimates$ATE

B = 500
estimates = data.frame(matrix(NA, nrow = B, ncol = 3))
for(b in 1:B){
  
  bootIndices = sample(1:n, replace = T)
  bootData = data[bootIndices,]
  X = subset(bootData, select = -c(studyid, any_violence, any_arrest))
  X = data.frame(model.matrix(~ . -1 , X))
  X1 = X0 = X
  X1$treatment = 1
  X0$treatment = 0
  newdata = rbind(X, X1, X0)
  Qinit = SuperLearner(Y=bootData$any_arrest, X=X, newX=newdata, SL.library=SL.library, family="binomial")
  QbarAW = Qinit$SL.predict[1:n]
  Qbar1W = Qinit$SL.predict[(n+1): (2*n)] 
  Qbar0W = Qinit$SL.predict[(2*n+1): (3*n)] 
  PsiHat.SS.b = mean(Qbar1W - Qbar0W)
  
  W = subset(X, select = -c(treatment))
  gHatSL = SuperLearner(Y=bootData$treatment, X=W, SL.library=SL.library, family="binomial")
  gHat1W = gHatSL$SL.predict
  gHat0W = 1- gHat1W
  
  PsiHat.IPTW.b = mean(bootData$treatment*bootData$any_arrest/gHat1W) - mean((1-bootData$treatment)*bootData$any_arrest/gHat0W)
  
  H.AW = (2*bootData$treatment-1)/ gHatAW
  H.1W = 1/gHat1W 
  H.0W = -1/gHat0W
  
  logitUpdate = glm(bootData$any_arrest ~ -1 +offset(qlogis(QbarAW)) + H.AW, family=binomial) 
  eps = logitUpdate$coef
  
  QbarAW.star = plogis(qlogis(QbarAW)+ eps*H.AW)
  Qbar1W.star = plogis(qlogis(Qbar1W)+ eps*H.1W)
  Qbar0W.star = plogis(qlogis(Qbar0W)+ eps*H.0W) 
  PsiHat.TMLE.b = mean(Qbar1W.star - Qbar0W.star) 
  estimates[b,] = c(PsiHat.SS.b, PsiHat.IPTW.b, PsiHat.TMLE.b)
}

colnames(estimates) = c("SimpSubps", "IPTW", "TMLE")

summary(estimates)
