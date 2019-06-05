library(monomvn)
IT<-10000
BI<-5000

ridgegibbs<-function(IT,BI){
  ridBayes<-bridge(matrix(c(day_before,day_2before,day_3before,day_4before,
                            day_5before,day_6before),ncol=6),day_now,T=IT,RJ=FALSE) #Bayesian ridge regression
  summary(ridBayes)
  
  est_beta1<-ridBayes$beta[,1]; est_beta2<-ridBayes$beta[,2]; est_beta3<-ridBayes$beta[,3]; est_beta4<-ridBayes$beta[,4];
  est_beta5<-ridBayes$beta[,5]; est_beta6<-ridBayes$beta[,6]; 
  
  est_int<-ridBayes$mu
}

par(mfrow=c(1,2))
plot(est_int[BI:IT], xlab = 'Index', ylab = ' Estimate for the intercept', main='Trace of Intercept',type='l')
abline(h=true_int, col='red')
hist(est_int[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Intercept')
abline(v=true_int, col='red')

par(mfrow=c(2,3))
hist(est_beta1[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Beta1')
abline(v=true_beta1, col='red')

hist(est_beta2[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Beta2')
abline(v=true_beta2, col='red')

hist(est_beta3[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Beta3')
abline(v=true_beta3, col='red')

hist(est_beta4[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Beta4')
abline(v=true_beta4, col='red')

hist(est_beta5[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Beta5')
abline(v=true_beta5, col='red')

hist(est_beta6[BI:IT],xlab = 'Index',ylab='Frequence', main='Histogram of the Beta6')
abline(v=true_beta6, col='red')


mayeight<-mean(est_int[BI:IT])*array(1,length(day_before))+mean(est_beta1[BI:IT])*day_before+mean(est_beta2[BI:IT])*day_2before+mean(est_beta3[BI:IT])*day_3before+mean(est_beta4[BI:IT])*day_4before+mean(est_beta5[BI:IT])*day_5before+mean(est_beta6[BI:IT])*day_6before

boxplot(mayeight)

n=length(day_before)
b=floor(n^(seq(0.01,0.7,by=0.01)))

#window estimator
w<-seq(from=1,to=0,by=-0.0009165903)
allasymvars<-array(0,n+1)
for (i in 1:n+1){
  asympvarest<-array(0,n+1-i)
  for(j in 1:n+1-i){
    asympvarest[j]<-(mayeight[i+j]-mean(na.omit(mayeight)))*(mayeight[i]-mean(na.omit(mayeight)))
  }
  if (i==1){
  allasymvars[i]<-w[i]/n*sum(na.omit(as.matrix(asympvarest))) 
  } else{
    allasymvars[i]<-2*w[i]/n*sum(na.omit(as.matrix(asympvarest))) 
  }

}
estimate<-sum(allasymvars)
estimate<-290*n
confidenceinterval<-cbind(mayeight-1.96*sqrt(abs(array(estimate,length(mayeight))))/sqrt(n),mayeight+sqrt(abs(array(estimate,length(mayeight))))/sqrt(n)*1.96)
truevalue<-as.numeric(day_now)
par(mfrow=c(1,1))

controle<-array(0,n)
for (i in 1:(n-1)){
  if (truevalue[i]>confidenceinterval[i,1] & truevalue[i]<confidenceinterval[i,2]){
    controle[i]=1
  }else{
    controle[i]=0
  }
}
sum(controle)/length(controle)

plot(truevalue)
lines(truevalue)
lines(confidenceinterval[,1][2:n],col='red')
lines(confidenceinterval[,2][2:n],col='green')

