rm(list=ls())

sigma = 1
phi=-0.5
n=10000;
bopti=floor((2*abs(phi)/(1-phi^2))^(2/3)*n^(1/3));
assympvar<-sigma^2/(1-phi)^2
experiment<-function(n,b,k){
  bmasvar <- function(n, b, k, sigma, phi) {
    X<-array(0,n)
    X[1]<-rnorm(1,0,sigma/(1-phi^2))
    
    Z<-rnorm(n-1,0,sigma)
    
    i=0
    for (i in c(2:n)){
      X[i]=phi*X[i-1]+Z[i-1];
    }
    
    average<-array(mean(X),k)
    
    batchmeans<-.colMeans(X,b,k)
    
    
    varb<-1/(k-1)*sum((batchmeans-average)^2)
    assympvarest<-varb*b
  }
  
  estimates<-array(0,length(b))
  for (i in 1:length(b)){
    estimates[i]<-bmasvar(n[i],b[i],k[i],sigma,phi)
  }
  estimates;
}

m = 100;
b=floor(n^seq(from= 0, to = 1, by = 0.01))
k = floor (n/b) 
n = b*k

totaal<-matrix(0,nrow=m,ncol=length(b))
means<-array(0,length(b))
variances<-array(0,length(b))
smallestbp<-array(0,length(b))

for (j in 1:m){
  totaal[j,]<-experiment(n,b,k)
}

for (i in 1:length(b)){
  means[i]<-mean(na.omit(totaal[,i]))
  variances[i]<-var(na.omit(totaal[,i]))
  smallestbp[i]<-(4-mean(totaal[,i]))^2+var(totaal[,i])
}

bestind<-which.min(smallestbp)
colnames(totaal)<-seq(from= 0, to = 1, by = 0.01)

bias<- 2*phi/((1-phi)*(1-phi^2)*b)
variance<-2/((1-phi)^4)*b/n

layout(matrix(c(1,1,1,1,1,1,2,2), 4, 2, byrow = TRUE))
par(mfrow=(c(1,1)))
boxplot(totaal,main = paste('b0 = ',as.character(bopti),', optimal alpha = ', as.character(seq(from= 0, to = 1, by = 0.01)[36]), 'phi = 0.9'),xlab='Alpha', ylab='Estimator for asymptotic variance')
abline(h=assympvar,col='red')
plot(seq(from= 0, to = 1, by = 0.01),variance,col='blue',type='l',xlab='Alpha',ylab='Size of the Bias/Variance',main='Size of the Bias and Variance for different Alpha')
lines(seq(from= 0, to = 1, by = 0.01),bias^2,type='l',col='red', xlab='Alpha',ylab='Size of the Bias^2/Variance')


plot(seq(from= 0, to = 1, by = 0.01)[1:70],smallestbp[1:70],main='Optimal MSE (LSE) for different alpha',ylab='MSE',xlab='Alpha')
abline(h=min(smallestbp[1:70]),col='red')
seq(from= 0, to = 1, by = 0.01)[which.min(smallestbp[1:70])]
