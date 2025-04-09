data<-read.csv("PdamRwarming.csv",header=T)
data
summary(data)
attach(data)
hist(total)
hist(log(total))
#defining vars and separating data
Dcol<-subset(data,sym=="D")
Dcol
Ccol<-subset(data,sym=="C")
Ccol
Dcolfeb<-subset(Dcol,date=="feb")
Dcolapr<-subset(Dcol,date=="apr")
Dcoljun<-subset(Dcol,date=="jun")
Ccolfeb<-subset(Ccol,date=="feb")
Ccolapr<-subset(Ccol,date=="apr")
Ccoljun<-subset(Ccol,date=="jun")

#data exploration, assumption testing
#testing if logtotals grouped by sym and month are normal
shapiro.test((log(Dcolfeb$total)))
shapiro.test((log(Dcolapr$total)))
shapiro.test((log(Dcoljun$total)))
shapiro.test((log(Ccolfeb$total)))
shapiro.test((log(Ccolapr$total)))
shapiro.test((log(Ccoljun$total)))
#logtotals grouped by sym and month are all normal!!
#histograms of logtotals to look for outliers
hist((log(Dcolfeb$total)))
hist((log(Dcolapr$total)))
hist((log(Dcoljun$total)))
hist((log(Ccolfeb$total)))
hist((log(Ccolapr$total)))
hist((log(Ccoljun$total)))
#boxplots of logtotals to look for outliers
boxplot((log(Dcolfeb$total)))
boxplot((log(Dcolapr$total))) #possible outlier colony 775, apr
(log(Dcolapr$total)-mean(log(Dcolapr$total)))/sd(log(Dcolapr$total))
#colony 775, apr--value is 2.09 standard deviations above mean, not removed from dataset
boxplot((log(Dcoljun$total)))
boxplot((log(Ccolfeb$total)))
boxplot((log(Ccolapr$total)))
boxplot((log(Ccoljun$total)))  #possible outlier colony 791, jun
(log(Ccoljun$total)-mean(log(Ccoljun$total)))/sd(log(Ccoljun$total))  #colony 791, jun--value is 2.92 standard deviations above mean
#DID NOT REMOVE OUTLIERS because no real justification

plot(log(total)~time)


#Changes over time
changesdata<-read.csv("PdamRwarming_changes.csv")
attach(changesdata)
changesdata
hist(relchange)
boxplot(febtotal) 
hist(log(febtotal))		#DOES X VARIABLE NEED TO BE NORMAL???
plot(relchange~log(febtotal))
fullmodel<-lm(relchange~log(febtotal)*sym)
summary(fullmodel) #interaction between initial density and symbiont n.s., so remove
ancovamodel<-lm(relchange~log(febtotal)+sym) #"complete ANCOVA model"
summary(ancovamodel)
plot(ancovamodel$residuals~ancovamodel$fitted.values)
#now fit two models separately to calculate adjusted means
Ccoldata<-changesdata[sym=="C",]
Dcoldata<-changesdata[sym=="D",]
Ccolmodel<-lm(relchange~log(febtotal),data=Ccoldata)
summary(Ccolmodel)
Dcolmodel<-lm(relchange~log(febtotal),data=Dcoldata)
summary(Dcolmodel)

#examine residuals of the two linear models
plot(residuals(Ccolmodel)~Ccolmodel$fitted.values)
plot(residuals(Ccolmodel)~log(febtotal),data=Ccoldata)
boxplot(residuals(Ccolmodel))
hist(residuals(Ccolmodel))
qqnorm(residuals(Ccolmodel))
qqline(residuals(Ccolmodel))
shapiro.test(residuals(Ccolmodel))
cdC<-cooks.distance(Ccolmodel)
plot(cdC~log(febtotal),data=Ccoldata,ylim=c(0,.3))
abline(h=4/21,lty=2)
#Ccolmodel residuals look good
plot(residuals(Dcolmodel)~Dcolmodel$fitted.values)
plot(residuals(Dcolmodel)~log(febtotal),data=Dcoldata)
boxplot(residuals(Dcolmodel))
hist(residuals(Dcolmodel))
qqnorm(residuals(Dcolmodel))
qqline(residuals(Dcolmodel))
shapiro.test(residuals(Dcolmodel))
cdD<-cooks.distance(Dcolmodel)
plot(cdD~log(febtotal),data=Dcoldata,ylim=c(0,.35))
abline(h=4/30,lty=2)
identify(log(febtotal)[sym=="D"],cdD,Dcoldata$colony)


#now calculate adjusted means using equations of each individual model with overall mean value for X (mean log(febtotal))

Ccolmodel$coefficients; Dcolmodel$coefficients
summary(log(febtotal))
Dcoldata
mean(relchange[sym=="C"]) #unadjusted mean
mean(relchange[sym=="D"]) #unadjusted mean
summary(relchange)
#plotting…
plot(relchange~log(febtotal),data=Ccoldata,pch=1,xlim=c(-4.1,-1.4),ylim=c(-0.5,1.2))
points(relchange~log(febtotal),data=Dcoldata,pch=19)
lines(predict(Ccolmodel)~log(febtotal),lty=4,data=Ccoldata)
lines(predict(Dcolmodel)~log(febtotal),data=Dcoldata)
abline(h=0,lty=2)



#model relating juntotal to febtotal
plot(juntotal~febtotal,data=changesdata) #looks like data needs transformation
plot(log(juntotal)~log(febtotal),data=changesdata)     #transformed data look appropriate for regression
summary(lm(log(juntotal)~log(febtotal)))
fullmodel<-lm(log(juntotal)~log(febtotal)*sym,data=changesdata)
summary(fullmodel)	#interaction n.s., go to complete ancova model
ancovamodel<-lm(juntotal~febtotal+sym,data=changesdata)
summary(ancovamodel)	#sym coefficient significant, log(febtotal) coeff (=slope) needs to be tested against 1…
						#t=(0.1580-1)/0.1523=-5.529 = significant (compared to t[.95,52]=1.675ish)
#make individual models for C and D
Cmodel<-lm(log(juntotal)~log(febtotal),data=changesdata[sym=="C",])
summary(Cmodel)	#t=(0.1580-1)/.1347=-6.251   SLOPE IS LESS THAN 1
Dmodel<-lm(log(juntotal)~log(febtotal),data=changesdata[sym=="D",])
summary(Dmodel)	#t=(0.4295-1)/.1252=-4.557	SLOPE IS LESS THAN 1
#examine C model residuals and leverage
plot(residuals(Cmodel)~Cmodel$fitted.values)
plot((residuals(Cmodel)-mean(residuals(Cmodel)))/sd(residuals(Cmodel))~Cmodel$fitted.values,ylim=c(-3,3),xlab="Fitted values",ylab="Standardized Residuals")
abline(h=0,lty=2)
qqnorm(residuals(Cmodel))
qqline(residuals(Cmodel))
hist(residuals(Cmodel))
shapiro.test(residuals(Cmodel))
plot(cooks.distance(Cmodel)~Cmodel$fitted.values)
abline(h=4/22,lty=2)
#examine D model residuals
plot(residuals(Dmodel)~Dmodel$fitted.values)
plot((residuals(Dmodel)-mean(residuals(Dmodel)))/sd(residuals(Dmodel))~Dmodel$fitted.values,ylim=c(-3,3),xlab="Fitted values",ylab="Standardized Residuals")
abline(h=0,lty=2)
qqnorm(residuals(Dmodel))
qqline(residuals(Dmodel))
hist(residuals(Dmodel))
shapiro.test(residuals(Dmodel))
plot(cooks.distance(Dmodel)~Dmodel$fitted.values)
abline(h=4/30,lty=2)
#plotting
plot(log(juntotal)~log(febtotal),data=changesdata[sym=="D",],xlim=c(-4,-1),ylim=c(-4,-1),pch=20,xlab="Feb. Density",ylab="June Density")
points(log(juntotal)~log(febtotal),data=changesdata[sym=="C",])
lines(predict(Dmodel)~log(febtotal),data=changesdata[sym=="D",])
lines(predict(Cmodel)~log(febtotal),data=changesdata[sym=="C",],lty=2)
abline(coef=c(0,1),lty=1)