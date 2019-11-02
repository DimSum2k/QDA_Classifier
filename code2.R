library(MASS)

#Question 1 -----------------------------------------------------

load("Xtest.RData") #chargement des donnees
load("Ytest.RData")
load("Xtrain.RData")
load("Ytrain.RData")
obs = length(Ytrain) + length(Ytest)
print(c("Number of observations: ", obs))
print(c("Number of observations in train set: ",length(Ytrain))) #numbers of images * 40 000
print(c("Number of observations in test set: ",length(Ytest)))
print(c("Dimension of each observation: ",length(Xtrain[1,])))
par(mfrow=c(2,2), mai=c(0.8,0.1,0.1,0.1)) #display images on a grid
for (i in c(5,22,45,47)) {
  mat = matrix(rev(Xtest[i,]), nrow=200, byrow=TRUE)#resize matrix, #dim(mat)=200*200
  string = "Image of a"
  if (Ytest[i]==1) string = paste(string,"cat, labelled",toString(1))
  else string = paste(string,"dog, labelled ",toString(0))
  
  image(mat,col = grey(seq(0, 1, length = 256)), xaxt= "n", yaxt= "n", xlab=string  ) #Display image
}

print(c("Percentage of positive targets (cats) in train set",sum(Ytrain==1)/length(Ytrain)))
print(c("Percentage of positive targets (cats) in test set",sum(Ytest==1)/length(Ytest)))


#Question 2 -----------------------------------------------------


compute_PCA <- function(Xtrain,Xtest) {
  X = rbind(Xtrain,Xtest) #concatenation
  X = scale(X, center = TRUE, scale = FALSE)  #centrage 
  PCA = svd(X,nu=0,nv=15) #nv=15 on ne conserve que 15 composantes
  C = X %*% PCA$v
  Ctrain = C[1:315,]
  Ctest = C[316:dim(C)[1],]
  save(Ctrain, file = "Ctrain.RData") #sauvegarde des donnees reduites
  save(Ctest, file = "Ctest.RData")
  
  return(sum(PCA$d[1:15]**2)/sum(PCA$d**2)) #variance expliquee
}

compute_PCA(Xtrain,Xtest)


#Question 9 -----------------------------------------------------

load("Ctrain.RData")
load("Ctest.RData")

computeML <- function(C, Y){
  n = length(Y)
  N1 = sum(Y==1)
  C0 = C[Y==0,]
  C1 = C[Y==1,]
  
  p_hat = N1/n
  mu_hat0 = colMeans(C0)
  mu_hat1 = colMeans(C1)
  C0_centered = sweep(C0,2,mu_hat0) #susbtract mu_hat0 to C0 
  C1_centered = sweep(C1,2,mu_hat1)
  sigma_hat0 = t(C0_centered)%*%C0_centered/(n-N1)
  sigma_hat1 = t(C1_centered)%*%C1_centered/N1
  
  out = list(p_hat,mu_hat0,mu_hat1,sigma_hat0,sigma_hat1)
  
  return(out)
}

ML = computeML(Ctrain,Ytrain) 

#QDA from MASS package
qda.model = qda(Ctrain,Ytrain)

ML = computeML(Ctrain,Ytrain) 
#Comparaisons
cat("Estimation of p : ", ML[[1]])
cat("QDA estimation of p : ", qda.model$prior[2])
cat("Estimation of log(det(sigma0)) : ", log(det(ML[[4]])))
cat("QDA estimation of log(det(sigma0)) : ", qda.model$ldet[1])
cat("Estimation of log(det(sigma1)) : ", log(det(ML[[5]])))
cat("QDA estimation of log(det(sigma1)) : ", qda.model$ldet[2])
cat("Estimation of mu0[1:4] : ",  ML[[2]][1:4])
cat("QDA estimation of mu0[1:4] : ",  qda.model$means[1,1:4])
cat("Estimation of mu1[1:4] : ",  ML[[3]][1:4])
cat("QDA estimation of mu1[1:4] : ",  qda.model$means[2,1:4])


#Question 14 -----------------------------------------------------


computeLogRatio <- function(cvect,p,mu0,mu1,Sigma0,Sigma1) {
  
  logratio = (0.5*(-log(det(Sigma1)) + log(det(Sigma0)) 
                   - t(cvect-mu1)%*%ginv(Sigma1)%*%(cvect-mu1)
                   + t(cvect-mu0)%*%ginv(Sigma0)%*%(cvect-mu0)) 
              + log(p) - log(1-p))
  
  return(logratio) 
}

computePred <- function(C,p,mu0,mu1,Sigma0,Sigma1) {
  
  toapply <- function(cvect,p,mu0,mu1,Sigma0,Sigma1) {
    return(as.integer((computeLogRatio(cvect,p,mu0,mu1,Sigma0,Sigma1)>0)))
  }
  
  pred = apply(C, MARGIN = 1, FUN = toapply, p = p, mu0 = mu0, 
               mu1 = mu1, Sigma0 = Sigma0, Sigma1 = Sigma1)
  
  return(pred)
}


#Question 15 -----------------------------------------------------

stats = computeML(Ctrain,Ytrain)
prediction = computePred(Ctest,stats[[1]],stats[[2]],stats[[3]],stats[[4]],stats[[5]])

prediction

#Erreur de prediction :
sum(prediction==Ytest)/length(Ytest)

computeQDA <- function(Ctrain, Ctest, Ytrain, Ytest){
  qda.model = qda(Ctrain,Ytrain)
  pred = predict(qda.model, Ctest)
  print(pred$class)
  return(sum(pred$class==Ytest)/length(Ytest))
}

computeQDA(Ctrain, Ctest, Ytrain, Ytest)


#Bonus -----------------------------------------------------


computeMLlin <- function(C, Y){
  #Calcul les statistiques associees a la LDA
  n = length(Y)
  N1 = sum(Y==1)
  C0 = C[Y==0,]
  C1 = C[Y==1,]
  p_hat = N1/n
  mu_hat0 = colMeans(C0)
  mu_hat1 = colMeans(C1)
  C0_centered = sweep(C0,2,mu_hat0)
  C1_centered = sweep(C1,2,mu_hat1)
  Sigma = (t(C0_centered)%*%C0_centered + t(C1_centered)%*%C1_centered)/(n-2)
  
  out = list(p_hat,mu_hat0,mu_hat1,Sigma)
  return(out)
}

computeLogRatiolin <- function(c,p,mu0,mu1, Sigma) {
  #Calcul les log ratio associees a la LDA 
  
  logratio = (log(p) - log(1-p) - 0.5*t(mu1 + mu0)%*%ginv(Sigma)%*%(mu1 - mu0) 
              + t(c)%*%ginv(Sigma)%*% (mu1-mu0))
  
  return(logratio)
}

computePredlin <- function(C,p,mu0,mu1,Sigma) {
  
  nbobs = dim(C)[1]
  pred = rep(0,nbobs)
  for (i in 1:nbobs) {
    
    logratio = computeLogRatiolin(C[i,],p,mu0,mu1,Sigma) 
    predi = logratio>0
    pred[i] = predi
  }
  
  return(pred) 
}

computeLDA <- function(Ctrain, Ctest, Ytrain, Ytest){
  lda.model = lda(Ctrain,Ytrain)
  pred = predict(lda.model, Ctest)
  print(pred$class)
  return(sum(pred$class==Ytest)/length(Ytest))
}


statsbis = computeMLlin(Ctrain,Ytrain)
prediction = computePredlin(Ctest,statsbis[[1]],statsbis[[2]],statsbis[[3]],statsbis[[4]])
prediction
sum(prediction==Ytest)/length(Ytest)
computeLDA(Ctrain, Ctest, Ytrain, Ytest)


