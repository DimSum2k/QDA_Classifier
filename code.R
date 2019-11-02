load("Xtest.RData")
load("Ytest.RData")
load("Xtrain.RData")
load("Ytrain.RData")
setwd("/Users/dimitri/Library/Mobile Documents/com~apple~CloudDocs/M1ENSAE/Semestre1/Statistique1/QDA_Classifier")
getwd()
list.files()

compute_PCA <-function() {
  #load("/Xtest.RData") #chargement des donnees
  #load("/Xtrain.RData")
  X = rbind(Xtrain,Xtest) #concatenation
  X = scale(X, center = TRUE, scale = FALSE)  #centrage 
  PCA = svd(X,nu=0,nv=15) #nv=15 on ne conserve que les 15 PCA
  C = X %*% PCA$v
  Ctrain = C[1:315,]
  Ctest = C[316:dim(C)[1],]
  save(Ctrain, file = "Ctrain.RData") #sauvegarde des donnees reduites
  save(Ctest, file = "Ctest.RData")
  
  return(sum(PCA$d[1:15]**2)/sum(PCA$d**2)) #variance expliquee
}

compute_PCA()
load("Ctrain.RData")
load("Ctest.RData")

computeML <- function(C, Y){
  n = length(Y)
  N1 = sum(Y==1)
  C1 = C[Y==1,]
  C0 = C[Y==0,]
  
  pi_hat = N1/n
  mu_hat1 = colMeans(C1)
  mu_hat0 = colMeans(C0)
  sigma_hat1 = t(sweep(C1,2,mu_hat1))%*%sweep(C1,2,mu_hat1)/N1
  sigma_hat0 = t(sweep(C0,2,mu_hat0))%*%sweep(C0,2,mu_hat0)/(n-N1)
  
  
  out = list(pi_hat,mu_hat0,mu_hat1,sigma_hat0,sigma_hat1)
  return(out)
}

computeLogRatio <- function(c,pi,mu0,mu1,Sigma0,Sigma1) {
  logratio = 0.5*(-log(det(Sigma1)) + log(det(Sigma0)) - t(matrix(c)-mu1)%*%Sigma1%*%(matrix(c)-mu1) + t(matrix(c)-mu0)%*%Sigma0%*%(matrix(c)-mu0)) + log(pi) - log(1-pi)

  return(logratio)
}


computePred <- function(Ctrain,Ctest,Ytrain,Ytest) {
  stat = computeML(Ctrain,Ytrain)
  score = 0
  for (i in 1:length(Ytest)) {
    logratio = computeLogRatio(Ctest[i,],stat[[1]],stat[[2]],stat[[3]],stat[[4]],stat[[5]])
    if (logratio>0) {result =1}
    else {result = 0}
    #print(c("Pred ",i, "is : ", result, "True result is : ", Ytest[i]))
    if (result==Ytest[i]) {score = score + 1}
    
  }
  
  return(score/length(Ytest))
}



computeMLlin <- function(C, Y){
  n = length(Y)
  N1 = sum(Y==1)
  C1 = C[Y==1,]
  C0 = C[Y==0,]
  
  pi_hat = N1/n
  mu_hat1 = colMeans(C1)
  mu_hat0 = colMeans(C0)
  Sigma = (t(sweep(C1,2,mu_hat1))%*%sweep(C1,2,mu_hat1) + t(sweep(C0,2,mu_hat0))%*%sweep(C0,2,mu_hat0))/(n-2)
  
  out = list(pi_hat,mu_hat0,mu_hat1,Sigma)
  return(out)
}
  
  computeLogRatiolin <- function(c,pi,mu0,mu1, Sigma) {
    logratio = log(pi/(1-pi)) - 0.5*t(mu1 + mu0)%*%Sigma%*%(mu1 - mu0) + t(matrix(c))%*%Sigma%*%(mu1-mu0)
    return(logratio)
  }


computePredlin <- function(Ctrain,Ctest,Ytrain,Ytest) {
  stat = computeMLlin(Ctrain,Ytrain)
  score = 0
  for (i in 1:length(Ytest)) {
    logratio = computeLogRatiolin(Ctest[i,],stat[[1]],stat[[2]],stat[[3]],stat[[4]])
    if (logratio>0) {result =1}
    else {result = 0}
    #print(c("Pred ",i, "is : ", result, "True result is : ", Ytest[i]))
    if (result==Ytest[i]) {score = score + 1}
    
  }
  
  return(score/length(Ytest))
}

computePred(Ctrain,Ctest,Ytrain,Ytest)
computePredlin(Ctrain,Ctest,Ytrain,Ytest)
