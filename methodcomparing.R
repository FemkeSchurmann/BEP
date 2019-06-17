# library(astsa)
  
phi = 0.1
sigma = 1
n = seq(from=200,to=10000,by=200)
b=n^0.5
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

boxplot(totaalbatch, main = 'phi=0.1 - Batch means method')
as.var.true = sigma^2/(1-phi)^2
abline(h=as.var.true,col='red')

boxplot(totaalblocks, main = 'phi=0.1 - Blocked means method')
as.var.true = sigma^2/(1-phi)^2
abline(h=as.var.true,col='red')

boxplot(t(totaal), main = 'phi=0.1 - Without window estimator')
abline(h=as.var.true,col='red')

boxplot(t(totaalwith), main = 'phi=0.1 - With window estimator')
abline(h=as.var.true,col='red')

