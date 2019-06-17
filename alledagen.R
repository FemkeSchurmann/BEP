rm(list=ls())
data1<-data.frame(`KNMI_20190501`)
library(monomvn)
data1$FHVEC[which(data1$FF>600,arr.ind = TRUE)]<-NA
rotterdam<-na.omit(data1)

day_now<-rotterdam$FHVEC[7:length(rotterdam$FHVEC)]
day_before<-rotterdam$FHVEC[6:(length(rotterdam$FHVEC)-1)]
day_2before<-rotterdam$FHVEC[5:(length(rotterdam$FHVEC)-2)]
day_3before<-rotterdam$FHVEC[4:(length(rotterdam$FHVEC)-3)]
day_4before<-rotterdam$FHVEC[3:(length(rotterdam$FHVEC)-4)]
day_5before<-rotterdam$FHVEC[2:(length(rotterdam$FHVEC)-5)]
day_6before<-rotterdam$FHVEC[1:(length(rotterdam$FHVEC)-6)]



regr_ln<-lm(day_now~day_before+day_2before+day_3before+day_4before+day_5before+day_6before)
summary(regr_ln)

anova(regr_ln)

par(mfrow=c(2,3))
plot(day_now,day_before, main = 'Impact of 1 day before on the windspeeds today', ylab='1 Day Before', xlab = 'Today')
plot(day_now,day_2before, main = 'Impact of 2 days before on the windspeeds today', ylab='2 Days Before', xlab = 'Today')
plot(day_now,day_3before, main = 'Impact of 3 days before on the windspeeds today', ylab='3 Days Before', xlab = 'Today')
plot(day_now,day_4before, main = 'Impact of 4 days before on the windspeeds today', ylab='4 Days Before', xlab = 'Today')
plot(day_now,day_5before, main = 'Impact of 5 days before on the windspeeds today', ylab='5 Days Before', xlab = 'Today')
plot(day_now,day_6before, main = 'Impact of 6 days before on the windspeeds today', ylab='6 Days Before', xlab = 'Today')

#ridge regression
rid<-regress(matrix(c(day_before,day_2before,day_3before,day_4before,
                      day_5before,day_6before),ncol=6), day_now, method='ridge')
summary(rid)
true_int<-rid$b[1]; true_beta1<-rid$b[2]; true_beta2<-rid$b[3]; true_beta3<-rid$b[4]
true_beta4<-rid$b[5]; true_beta5<-rid$b[6]; true_beta6<-rid$b[7]


#Gibbssampler ridge
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



#asymptoticcovariancematrix
findasympvar<-function(X){
  n = length(X)
  b = n^0.3 #here we choose the value for alpha
  k = floor(n/b)
  n = b*k
  resized_X<-X[1:n]
  
  batchmeans<-.colMeans(resized_X,b,k)
  varb<-1/(k-1)*sum((mean(resized_X)-batchmeans)^2)
  assympvarest<-b*varb
}

asympvarvector<-array(1,7)

asympvarvector[1]<-findasympvar(est_int[BI:IT]) #asymptotic variance of the intercept
asympvarvector[2]<-findasympvar(est_beta1[BI:IT])
asympvarvector[3]<-findasympvar(est_beta2[BI:IT])
asympvarvector[4]<-findasympvar(est_beta3[BI:IT])
asympvarvector[5]<-findasympvar(est_beta4[BI:IT])
asympvarvector[6]<-findasympvar(est_beta5[BI:IT])
asympvarvector[7]<-findasympvar(est_beta6[BI:IT])

asympcovariancematrix<-sqrt(asympvarvector)%*%t(sqrt(asympvarvector))/6

#predict more days with function
X<-cbind(array(1,length(day_before)),day_before,day_2before,day_3before,day_4before,day_5before,day_6before)
beta<-rbind(mean(est_int),mean(est_beta1),mean(est_beta2),mean(est_beta3),mean(est_beta4),mean(est_beta5),mean(est_beta6))

estimates<-X%*%beta
secondmoment<-mean(estimates^2)
firstmoment<-mean(estimates)
firstmomentsquared<-mean(estimates)^2
sigma<-17.22
g<-sqrt(sigma^2+secondmoment-firstmomentsquared)
gwithout<-sqrt(secondmoment-firstmomentsquared)
asympvariancefirstmoment<-findasympvar(estimates)
asympvariancesecondmoment<-findasympvar(estimates^2)

fout<-g+sqrt(asympvariancefirstmoment/length(estimates))*(1-firstmoment/(2*g))+1/(2*g)*sqrt(asympvariancesecondmoment/length(estimates))
foutwithout<- gwithout+sqrt(asympvariancefirstmoment/length(estimates))*(1-firstmoment/(2*gwithout))+1/(2*gwithout)*sqrt(asympvariancesecondmoment/length(estimates))

findconfidenceinterval<-function(day){
  X<-cbind(1,day_before[day],day_2before[day],day_3before[day],day_4before[day],day_5before[day],day_6before[day])
  beta<-rbind(mean(est_int),mean(est_beta1),mean(est_beta2),mean(est_beta3),mean(est_beta4),mean(est_beta5),mean(est_beta6))
  estimate_day<-X%*%beta
  
 
  confinterval<-cbind(estimate_day-1.96*fout,estimate_day-1.96*foutwithout,estimate_day,estimate_day+1.96*foutwithout,estimate_day+1.96*fout)
  return(confinterval)
  }


allintervals<-matrix(0,ncol=5,nrow=length(day_now))
for (i in 1:length(day_now)){
  allintervals[i,]<-findconfidenceinterval(i)
}




par(mfrow=c(1,1))
plot(day_now)
lines(day_now)

lines(allintervals[,1],col='red')

lines(allintervals[,3],col='red')
lines(allintervals[,2],col='green')

times<-seq(as.Date('2016-05-01'),by='days',length=length(estimates))

plot(times, day_now, type = "l",xlab='Date', ylab='Average windspeed on that day', main='Confidence interval for the estimates', ylim=c(-10,120))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(times,rev(times)),c(allintervals[,1],rev(allintervals[,5])),col = "cornsilk3", border = FALSE)
polygon(c(times,rev(times)),c(allintervals[,2],rev(allintervals[,4])),col = "cornsilk2", border = FALSE)
lines(times, day_now, lwd = 0.1,col='black')

controle<-array(0,length(day_now))
for (i in 1:(length(day_now)-1)){
  if (day_now[i]>allintervals[i,2] & day_now[i]<allintervals[i,4]){
    controle[i]=1
  }else{
    controle[i]=0
  }
}
sum(controle)/length(controle)
