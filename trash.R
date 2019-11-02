#Trash for the project

computeLogRatio2 <- function(c,pi,mu0,mu1,Sigma0,Sigma1) {
  
  logratio = 0.5*(-log(det(Sigma1)) + log(det(Sigma0)) 
                  - t(matrix(c)-mu1)%*%Sigma1%*%(matrix(c)-mu1) 
                  + t(matrix(c)-mu0)%*%Sigma0%*%(matrix(c)-mu0)) 
  + log(pi) - log(1-pi)
  
  return(logratio)
}



computePred2 <- function(Ctrain,Ctest,Ytrain,Ytest) {
  stat = computeML(Ctrain,Ytrain)
  result = rep(0, length(Ytest))
  
  for (i in 1:length(Ytest)) {
    logratio = computeLogRatio(Ctest[i,],stat[[1]],stat[[2]],stat[[3]],stat[[4]],stat[[5]])
    result[i] = (logratio>0)
    #print(c("Pred ",i, "is : ", result, "True result is : ", Ytest[i]))
  }
  print(result)
  return(sum(result==Ytest)/length(Ytest))
}