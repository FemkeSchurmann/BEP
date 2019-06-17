# library(astsa)
rm(list=ls())  

phi = -0.5
sigma = 1
n = seq(from=200,to=10000,by=200)
b=n^(0.33)
m = 100

totaal<-matrix(0,ncol=m,nrow=length(n))
totaalwith<-matrix(0,ncol=m,nrow=length(n))
par(mfrow=c(1,1))
for (j in 1:m){
asvars<-array(0,length(n))
asvarswith<-array(0,length(n))
for (i in 1:length(n)){
w = sigma * rnorm(n[i],0,1)
x<-filter(w, filter = c(phi), method="recursive");
gamma<-as.matrix(acf(x,type="covariance")$acf)


#Window
asvars[i] = gamma[1] + 2 * sum(gamma[2:length(gamma)])
asvarswith[i]= gamma[1] + sum((1+cos(pi*seq(from=length(gamma),to=2,by=-1)*4/n[i]))*gamma[2:length(gamma)])

}
totaal[,j]<-asvars
totaalwith[,j]<-asvarswith
}

#batchmeans
bopti=floor((2*abs(phi)/(1-phi^2))^(2/3)*n^(1/3));

bias<- 2*phi/((1-phi)*(1-phi^2)*bopti)


experimentbatch<-function(n,b,k){
  bmasvar <- function(n, b, k, sigma, phi) {
  w = sigma * rnorm(n,0,1)
  x<-filter(w, filter = c(phi), method="recursive");
  
    average<-array(mean(x),k)
    
    batchmeans<-.colMeans(x,b,k)
    
    
    varb<-1/(k-1)*sum((batchmeans-average)^2)
    assympvarest<-varb*b
  
  }
  estimates<-array(0,length(b))
for (i in 1:length(b)){
    estimates[i]<-bmasvar(n[i],b[i],k[i],sigma,phi)
  }
  estimates;
}


k = floor (n/b) 
n = b*k

totaalbatch<-matrix(0,nrow=m,ncol=length(b))

for (j in 1:m){
  totaalbatch[j,]<-experimentbatch(n,b,k)
}

#blockedmeans

experiment<-function(n,b){
  bmasvar <- function(n, b, sigma, phi) {
    w = sigma * rnorm(n,0,1)
    x<-filter(w, filter = c(phi), method="recursive");
    
    average<-array(mean(x),n-b+1)
    
    blockedmeans<-array(0,n-b+1)
    for (j in 1:(n-b+1)){
      blockedmeans[j]<-mean(x[j:(j+b-1)]) #mean of each block
    }
    
    
    varb<-1/(n-b+1)*sum((blockedmeans-average)^2)
    assympvarest<-varb*b
  }
  
  estimates<-array(0,length(b))
  for (i in 1:length(b)){
    estimates[i]<-bmasvar(n[i],b[i],sigma,phi)
  }
  estimates;
}

totaalblocks<-matrix(0,nrow=m,ncol=length(b))

for (j in 1:m){
  totaalblocks[j,]<-experiment(n,b)
}

n = seq(from=200,to=10000,by=200)
rownames(totaal)<-n
rownames(totaalwith)<-n
colnames(totaalblocks)<-n
colnames(totaalbatch)<-n

par(mfrow=c(2,2))

as.var.true = sigma^2/(1-phi)^2
batchMSE<-mean((colMeans(totaalbatch)-as.var.true)^2)
boxplot(totaalbatch, main = paste('phi=-0.5 - Batch means method,','MSE = ', round(batchMSE,digits=4) ),ylim=c(0,1))
abline(h=as.var.true,col='red')

blockMSE<-mean((colMeans(totaalblocks)-as.var.true)^2)
boxplot(totaalblocks, main = paste('phi=-0.5 - Blocked means method,','MSE = ', round(blockMSE,digits=4) ),ylim=c(0,1))
as.var.true = sigma^2/(1-phi)^2
abline(h=as.var.true,col='red')

windowMSE<-mean((colMeans(totaal)-as.var.true)^2)
boxplot(t(totaal), main = paste('phi=-0.5 - Without window estimator,','MSE = ', round(windowMSE,digits=4) ),ylim=c(0,1))
abline(h=as.var.true,col='red')

windowwithMSE<-mean((colMeans(totaalwith)-as.var.true)^2)
boxplot(t(totaalwith), main = paste('phi=-0.5 - With window estimator,','MSE = ', round(windowwithMSE,digits=4) ),ylim=c(0,1))
abline(h=as.var.true,col='red')

