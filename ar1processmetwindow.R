rm(list = ls())
sigma = 1
phi=0.5
n=seq(from=100,to=10000, by=100);
bopti=floor((2*abs(phi)/(1-phi^2))^(2/3)*n^(1/3));



allestimates<-function(n){
  X<-array(0,n)
  X[1]<-rnorm(1,0,sigma/(1-phi^2))
  
  Z<-rnorm(n-1,0,sigma)
  
  for (i in c(2:n)){
    X[i]=phi*X[i-1]+Z[i-1];
  }
  
  average<-mean(X)
  alllagedautocovs<-array(0,n)
  
  for (j in 1:n){ #here we caculate all autocovariance lags except for the lag 0
    h=j  
    
    lagautocov<-array(0,n-h)
    for (i in 1:(n-h)){
      lagautocov[i]<-(X[i]-average)*(X[i+h]-average)
    }
    alllagedautocovs[j]<-2*(1-sqrt((h/n)))/n*sum(lagautocov)
  }
  
  
  
  lagautocov<-array(0,n)
  for (i in 1:(n)){
    lagautocov[i]<-(X[i]-average)*(X[i]-average)
  }
  zerocase<-mean(lagautocov)
  
  
  estimate_asympvar<-sum(na.omit(as.numeric(alllagedautocovs)))+zerocase
  return(estimate_asympvar)
}


  tijdelijk<-array(1,length(n))
  for (l in 1:length(n)){
    tijdelijk[l]<-allestimates(n[l])
  }



par(mfrow=c(1,1))
plot(n,tijdelijk, main='Estimates for the asymptotic variance for different n', ylab='Estimate', xlab='Size of n')
lines(n,tijdelijk)
 abline(h=4,col='red')
 