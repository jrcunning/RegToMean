bleaching<-read.csv("PdamRbleaching.csv")
bleaching
attach(bleaching)
hist(juntotal)
boxplot(juntotal)
shapiro.test(juntotal)
hist(log(juntotal))
shapiro.test(log(juntotal))
plot(relchange~log(juntotal))
plot(pchange~log(juntotal))
hist(Cbleach$relchange)
hist(Dbleach$relchange)
Cbleach<-subset(bleaching,sym=="C")
Dbleach<-subset(bleaching,sym=="D")
plot(relchange~log(juntotal),data=Cbleach,xlim=c(-3.3,-0.9),ylim=c(-1.5,.5))
points(relchange~log(juntotal),data=Dbleach,pch=20)
plot(relchange~log(juntotal),data=bleaching,pch=as.character(sym),)

fullmodel<-lm(relchange~log(juntotal)*sym)
summary(fullmodel)
ancovamodel<-lm(relchange~log(juntotal)+sym)
summary(ancovamodel)
#fit individual models
Cmodel<-lm(relchange~log(juntotal),data=Cbleach)
summary(Cmodel)
Dmodel<-lm(relchange~log(juntotal),data=Dbleach)
summary(Dmodel)
#examine C model residuals
plot(residuals(Cmodel)~Cmodel$fitted.values)
plot(residuals(Cmodel)~log(juntotal),data=Cbleach)
boxplot(residuals(Cmodel))
hist(residuals(Cmodel))
qqnorm(residuals(Cmodel))
qqline(residuals(Cmodel))
shapiro.test(residuals(Cmodel))
cdC<-cooks.distance(Cmodel)
summary(Cbleach)
plot(cdC~log(juntotal),data=Cbleach,ylim=c(0,.3))
abline(h=4/19,lty=2)
identify(log(juntotal)[sym=="C"],cdC,Cbleach$colony)
#colony 981 has D>4/n
#examine D model residuals
plot(residuals(Dmodel)~Dmodel$fitted.values)
plot(residuals(Dmodel)~log(juntotal),data=Dbleach)
boxplot(residuals(Dmodel))
hist(residuals(Dmodel))
qqnorm(residuals(Dmodel))
qqline(residuals(Dmodel))
shapiro.test(residuals(Dmodel))
cdD<-cooks.distance(Dmodel)
summary(Dbleach)
plot(cdD~log(juntotal),data=Dbleach,ylim=c(0,.35))
abline(h=4/23,lty=2)
identify(log(juntotal)[sym=="D"],cdD,Dbleach$colony)
#colony 796 has D>4/n


#plotting
summary(relchange)
summary(log(juntotal))
plot(relchange~log(juntotal,10),data=Cbleach,pch=1,xlim=c(-1.4331718,-0.390865),ylim=c(-1.4,0.3),xlab="log Pre-bleaching density",ylab="Relative symbiont change",cex=0.7)
points(relchange~log(juntotal,10),data=Dbleach,pch=19,cex=0.7)
lines(predict(Cmodel)~log(juntotal),data=Cbleach)
lines(predict(Dmodel)~log(juntotal),data=Dbleach)
abline(h=0,lty=2)
identify(log(juntotal),relchange,colony)

#calculating adjusted means
Cmodel$coefficients; Dmodel$coefficients
mean(log(juntotal))
mean(Cbleach$relchange) #unadjusted mean
mean(Dbleach$relchange) #unadjusted mean
Cmodel$coefficients[1]+Cmodel$coefficients[2]*mean(log(juntotal))
Dmodel$coefficients[1]+Dmodel$coefficients[2]*mean(log(juntotal))

#two-way model
model<-lm(relchange~propD*log(juntotal),data=data1)
summary(model)
reduced<-lm(relchange~propD+log(juntotal))
summary(reduced)
single<-lm(relchange~log(juntotal))
summary(single)
single2<-lm(relchange~propD)
summary(single2)
residuals(reduced)
hist(reduced$residuals)
plot(reduced$residuals~log(juntotal))
plot(reduced$residuals~predict(reduced))

##can i make a two-way graph of this?
plot(relchange~propD*log(juntotal))


#ANCOVA with jun as predictor and aug as response, by clade
plot(log(augtotal)~log(juntotal),data=Cbleach,xlim=c(-3.3,-0.8),ylim=c(-5.2,-1.2),xlab="log(Pre-bleaching Density)",ylab="log(Post-bleaching Density)")
points(log(augtotal)~log(juntotal),data=Dbleach,pch=20)
lines(predict(Cmodel)~log(juntotal),data=bleaching[sym=="C",])
lines(predict(Dmodel)~log(juntotal),data=bleaching[sym=="D",])
abline(coef=c(0,1))
ancova<-lm(log(augtotal)~log(juntotal)+sym)
summary(ancova)
#need to test if slope is less than 1!!!!! HOW???

#individual model
Dmodel<-lm(log(augtotal)~log(juntotal),data=bleaching[sym=="D",])
Cmodel<-lm(log(augtotal)~log(juntotal),data=bleaching[sym=="C",])
summary(model)
plot(log(augtotal)~log(juntotal),xlim=c(-3,-1),ylim=c(-5.5,-1),data=bleaching[sym=="D",])
lines(predict(model)~log(juntotal),data=bleaching[sym=="D",])
abline(coef=c(0,1))

plot(log(augtotal)~log(juntotal),xlim=c(-3.3,-1.5),ylim=c(-5.5,-1),data=bleaching[sym=="C",])
lines(predict(model)~log(juntotal),data=bleaching[sym=="C",])
abline(coef=c(0,1))

summary(log(juntotal))

summary(relchange[sym=="D"])