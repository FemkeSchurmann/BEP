n = 100;

sequence<-seq(from=-5,to=5,by=(10/n))


X<-runif(length(sequence),0,1)

layout(matrix(c(1,1,1,1,1,1,2,2), 4, 2, byrow = TRUE))
plot(sequence, dnorm(sequence,0,1), type='l', main='Standard Normal Distribution', ylab='Frequence',xlab=' ')
points(X,dnorm(X,0,1),col='red', lwd= '0.01',cex=0.01)

plot(X,array(1,n+1),cex=0.01,yaxt='n', col = 'red',main='Random generated variables', xlab = 'Random samples from uniform distribution',ylab=' ',xlim = c(-5,5))


n= seq(10, 100000, by = 1000)
m=1000

averages<- matrix(0,ncol=length(n),nrow=m)

for (i in 1:length(n)){
  for (j in 1:m){
    averages[j,i]<-mean((1/(2*pi)^(1/2))*exp(-runif(n[i],0,1)^2/2))
  }
}

colnames(averages)<-n

par(mfrow=c(1,1))
boxplot(averages, main='Estimates for the Expectation',xlab='Sample size',ylab='Value for the estimate')
abline(h=0.3413,col='red')

hist(est_int[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Intercept')
abline(v=true_int, col='red')
