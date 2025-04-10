library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)
library(splines)

#first lets simulate some data over fake time-steps - mean changes, but variance stays the same

###################################
#case one, increasing growth, but static growth rate (mean diff of 1 between each time step)

set.seed(42) 

par(mfrow=c(2,4))
set.seed(42) #to ensure repeatability; skip to vary results slightly each time
tA<-rnorm(1000,mean=10,sd=2)
hist(tA, xlab="size",main="Time A")
text(x = -1, y = 260, labels = "(A)", xpd = NA)
abline(v = mean(tA), col="red", lwd=3, lty=2)

set.seed(1) #to ensure repeatability; skip to vary results slightly each time
tB<-rnorm(1000,mean=11,sd=2)
hist(tB, xlab="size",main="Time B")
text(x = -1, y = 260, labels = "(B)", xpd = NA)
abline(v = mean(tB), col="red", lwd=3, lty=2)

set.seed(16) #to ensure repeatability; skip to vary results slightly each time
tC<-rnorm(1000,mean=12,sd=2)
hist(tC, xlab="size",main="Time C")
text(x = -1, y = 260, labels = "(C)", xpd = NA)
abline(v = mean(tC), col="red", lwd=3, lty=2)

set.seed(89) #to ensure repeatability; skip to vary results slightly each time
tD<-rnorm(1000,mean=13,sd=2)
hist(tD, xlab="size",main="Time D")
text(x = -1, y = 260, labels = "(D)", xpd = NA)
abline(v = mean(tD), col="red", lwd=3, lty=2)

plot(tB~tA, xlab = "size at Time A", ylab = "size at Time B")
text(x = -0.5, y = 21, labels = "(E)", xpd = NA)
summary(lm(tB~tA)) #UNCORRELATED
plot(tC~tB, xlab = "size at Time B", ylab = "size at Time C")
text(x = 0, y = 21.5, labels = "(F)", xpd = NA)
summary(lm(tC~tB)) #UNCORRELATED
plot(tD~tC, xlab = "size at Time C", ylab = "size at Time D")
text(x = 0, y = 22, labels = "(G)", xpd = NA)
summary(lm(tC~tD)) #UNCORRELATED
plot(tD~tA, xlab = "size at Time A", ylab = "size at Time D")
text(x = 0, y = 22, labels = "(H)", xpd = NA)
summary(lm(tD~tA)) #UNCORRELATED

grAstat<-tB-tA
hist(grAstat)
mean(grAstat) #close to 1 as desired

grAstat2<-tB/tA
hist(grAstat2)
mean(grAstat2)

grBstat<-tC-tB
hist(grBstat)
mean(grBstat) #close to 1 as desired

grBstat2<-tC/tB
hist(grBstat2)
mean(grBstat2) 

grCstat<-tD-tC
hist(grCstat)
mean(grCstat) #close to 1 as desired

grCstat2<-tD/tC
hist(grCstat2)
mean(grCstat2)

grDstat<-tD-tA
hist(grDstat)
mean(grDstat) #for a total of almost 3

grDstat2<-tD/tA
hist(grDstat2)
mean(grDstat2)

##now combine into dataframe
A<-as.data.frame(cbind("size"=tA,"growthRate"=grAstat,"growthRate2"=grAstat2,"InitialSize"=tA))
A$time<-rep("A",1000)
A$diff<-rep("A to A",1000)
head(A)

B<-as.data.frame(cbind("size"=tB,"growthRate"=grAstat, "growthRate2"=grAstat2,"InitialSize"=tA))
B$time<-rep("B",1000)
B$diff<-rep("A to B",1000)
head(B)

C<-as.data.frame(cbind("size"=tC,"growthRate"=grBstat, "growthRate2"=grBstat2,"InitialSize"=tB))
C$time<-rep("C",1000)
C$diff<-rep("B to C",1000)
head(C)

D<-as.data.frame(cbind("size"=tD,"growthRate"=grCstat, "growthRate2"=grCstat2,"InitialSize"=tC))
D$time<-rep("D",1000)
D$diff<-rep("C to D",1000)
head(D)

# Total<-as.data.frame(cbind("size"=tD,"growthRate"=grDstat,"InitialSize"=tA))
# Total$time<-rep("Total",1000)
# Total$diff<-rep("Total",1000)
# head(Total)


simStatic<-as.data.frame(rbind(A,B,C,D))
simStatic$stdGR<-(simStatic$growthRate/simStatic$InitialSize) #calculate standardized growth rate
simStatic$time<-factor(simStatic$time,ordered = T,levels=c("A","B","C","D")) #making time ordinal

#and generate plots and statistics

#For size over time
StatSize<-ggplot(simStatic, aes(x=time, y=size)) + geom_boxplot() +theme_bw()+ylab("Size")+xlab(NULL)+ggtitle("Static Growth Rate")+ geom_signif(comparisons = list(c("A", "B"),c("B", "C"),c("C", "D")), map_signif_level=TRUE) #simulated size data gets bigger over time steps on average, but no correlation between individuals across time-points
Mss<-aov(size~time,simStatic)
summary(Mss)
TukeyHSD(Mss,which="time") #size at each time-step sig diff from all others

#For growth over time
subStat1=subset(simStatic,time!="A") #3 pairwise differences for 4 time-points
StatGrowth<-ggplot(subStat1, aes(x=diff, y=growthRate)) + geom_boxplot() +theme_bw() +ylab("Growth") +xlab(NULL) 
Msgr<-aov(growthRate~diff,subStat1) 
summary(Msgr) #no change over time


#For standardized growth over time
StatGrowthStd<-ggplot(subStat1, aes(x=diff, y=stdGR)) + geom_boxplot() +theme_bw() +ylab("Std Growth") +xlab(NULL) + geom_signif(comparisons = list(c("A to B", "C to D")), map_signif_level=TRUE) #this is spurious significance based on table below - need to correct in final plot
Msgrs<-aov(stdGR~diff,subStat1) 
summary(Msgrs) #decreasing over time
TukeyHSD(Msgrs,which="diff") #sig difference bewteen first and last time windows


#for growth rate vs initial size v1
StatGrowthVSsize<-ggplot(subStat1, aes(x=InitialSize, y=growthRate)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+theme_bw() +ylab("Growth") +xlab(NULL) +ggtitle("Static Growth Rate")
summary(lm(growthRate~InitialSize,subStat1)) #significant negative relationship
par(mfrow=c(2,2))
plot(lm(growthRate~InitialSize,subStat1)) #satisfies assumptions of linear model
# see https://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/

StatGrowth2VSsize<-ggplot(subStat1, aes(x=InitialSize, y=growthRate2)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+theme_bw() +ylab("Growth") +xlab(NULL) +ggtitle("Static Growth Rate")
summary(lm(growthRate~InitialSize,subStat1)) #significant negative relationship
par(mfrow=c(2,2))
plot(lm(growthRate~InitialSize,subStat1)) #satisfies assumptions of linear model
# see https://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/



###now, plot standardized growth rates and see how relationship looks

StatStdGrowthVSsize<-ggplot(subStat1, aes(x=InitialSize, y=stdGR)) + geom_point(alpha=0.3) +geom_smooth()+theme_bw()+ylab("Std. Growth") +xlab(NULL) 
#strong negative non-linear?
par(mfrow=c(2,2))
plot(lm(stdGR~InitialSize,subStat1)) #yep, violates assumptions 


############################
#case two, increasing growth, with increasing growth rate (mean diff increases each time step)

par(mfrow=c(2,4))
set.seed(42) #to ensure repeatability; skip to vary results slightly each time
tA<-rnorm(1000,mean=10,sd=2)
hist(tA, xlab="size",main="Time A")
text(x = -1, y = 260, labels = "(A)", xpd = NA)
abline(v = mean(tA), col="red", lwd=3, lty=2)

tB<-rnorm(1000,mean=11,sd=2)
set.seed(1) #to ensure repeatability; skip to vary results slightly each time
hist(tB, xlab="size",main="Time B")
text(x = -3, y = 250, labels = "(B)", xpd = NA)
abline(v = mean(tB), col="red", lwd=3, lty=2)

tC<-rnorm(1000,mean=13,sd=2)
set.seed(16) #to ensure repeatability; skip to vary results slightly each time
hist(tC, xlab="size",main="Time C")
text(x = -0.5, y = 260, labels = "(C)", xpd = NA)
abline(v = mean(tC), col="red", lwd=3, lty=2)

tD<-rnorm(1000,mean=16,sd=2)
set.seed(89) #to ensure repeatability; skip to vary results slightly each time
hist(tD, xlab="size",main="Time D")
text(x = 3, y = 260, labels = "(D)", xpd = NA)
abline(v = mean(tD), col="red", lwd=3, lty=2)

plot(tB~tA, xlab = "size at Time A", ylab = "size at Time B")
text(x = -0.5, y = 21, labels = "(E)", xpd = NA)
summary(lm(tB~tA)) #UNCORRELATED
plot(tC~tB, xlab = "size at Time B", ylab = "size at Time C")
text(x = 1, y = 23, labels = "(F)", xpd = NA)
summary(lm(tC~tB)) #UNCORRELATED
plot(tD~tC, xlab = "size at Time C", ylab = "size at Time D")
text(x = 2, y = 25, labels = "(G)", xpd = NA)
summary(lm(tC~tD)) #UNCORRELATED
plot(tD~tA, xlab = "size at Time A", ylab = "size at Time D")
text(x = -2, y = 25, labels = "(H)", xpd = NA)
summary(lm(tD~tA)) #UNCORRELATED

grAinc<-tB-tA
hist(grAinc)
mean(grAinc) #nearly 1

grBinc<-tC-tB
hist(grBinc)
mean(grBinc) #nearly 2

grCinc<-tD-tC
hist(grCinc)
mean(grCinc) #nearly 3

grDinc<-tD-tA
hist(grDinc)
mean(grDinc) #for a total of 6

##now combine into dataframe
A<-as.data.frame(cbind("size"=tA,"growthRate"=grAinc,"InitialSize"=tA))
A$time<-rep("A",1000)
A$diff<-rep("A to A",1000)
head(A)

B<-as.data.frame(cbind("size"=tB,"growthRate"=grAinc,"InitialSize"=tA))
B$time<-rep("B",1000)
B$diff<-rep("A to B",1000)
head(B)

C<-as.data.frame(cbind("size"=tC,"growthRate"=grBinc,"InitialSize"=tB))
C$time<-rep("C",1000)
C$diff<-rep("B to C",1000)
head(C)

D<-as.data.frame(cbind("size"=tD,"growthRate"=grCinc,"InitialSize"=tC))
D$time<-rep("D",1000)
D$diff<-rep("C to D",1000)
head(D)


simInc<-as.data.frame(rbind(A,B,C,D))
simInc$stdGR<-(simInc$growthRate/simInc$InitialSize) #calculate standardized growth rate
simInc$time<-factor(simInc$time,ordered = T,levels=c("A","B","C","D")) #making time ordinal

#and plot
IncSize<-ggplot(simInc, aes(x=time, y=size)) + geom_boxplot() +theme_bw()+ylab(NULL)+xlab(NULL) +ggtitle("Incr. Growth Rate")+ geom_signif(comparisons = list(c("A", "B"),c("B", "C"),c("C", "D")), map_signif_level=TRUE)#simulated size data gets bigger over time steps - pseudo exponential
Mis<-aov(size~time,simInc)
summary(Mis)
TukeyHSD(Mis,which="time") #size at each time-step sig diff from all others

#For growth over time
subInc1=subset(simInc,time!="A") #3 pairwise differences for 4 time-points
IncGrowth<-ggplot(subInc1, aes(x=diff, y=growthRate)) + geom_boxplot() +theme_bw() +ylab(NULL) +xlab(NULL) + geom_signif(comparisons = list(c("A to B", "B to C"),c("B to C", "C to D")), map_signif_level=TRUE)
Migr<-aov(growthRate~diff,subInc1) 
summary(Migr) #sig diff over time
TukeyHSD(Migr,which="diff") #increasing (diff always positive)

#For standardized growth over time
IncGrowthStd<-ggplot(subInc1, aes(x=diff, y=stdGR)) + geom_boxplot() +theme_bw() +ylab(NULL) +xlab(NULL) + geom_signif(comparisons = list(c("A to B", "B to C"),c("B to C", "C to D")), map_signif_level=TRUE)
Migrs<-aov(stdGR~diff,subInc1) 
summary(Migrs) #decreasing over time
TukeyHSD(Migrs,which="diff") #sig difference bewteen first and last time windows

#for growth rate vs initial size 
IncGrowthVSsize<-ggplot(subInc1, aes(x=InitialSize, y=growthRate)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+theme_bw() +ylab(NULL) +xlab(NULL) +ggtitle("Increasing Growth Rate")
summary(lm(growthRate~InitialSize,subInc1)) #significant negative relationship


###now, plot standardized growth rates and see how relationship looks

IncStdGrowthVSsize<-ggplot(subInc1, aes(x=InitialSize, y=stdGR)) + geom_point(alpha=0.3) +geom_smooth()+theme_bw()+ylab(NULL) +xlab(NULL)
 #very strong negative (and nonlinear) relationship
par(mfrow=c(2,2))
plot(lm(stdGR~InitialSize,subInc1)) #yep, violates assumptions 


######################################################################################
#case three, increasing growth, with decreasing growth rate (mean diff decreases each time step)
par(mfrow=c(2,4))
set.seed(42) #to ensure repeatability; skip to vary results slightly each time
tA<-rnorm(1000,mean=10,sd=2)
hist(tA, xlab="size",main="Time A")
text(x = -1, y = 240, labels = "(A)", xpd = NA)
abline(v = mean(tA), col="red", lwd=3, lty=2)

set.seed(1) #to ensure repeatability; skip to vary results slightly each time
tB<-rnorm(1000,mean=13,sd=2)
hist(tB, xlab="size",main="Time B")
text(x = 2, y = 250, labels = "(B)", xpd = NA)
abline(v = mean(tB), col="red", lwd=3, lty=2)

set.seed(16) #to ensure repeatability; skip to vary results slightly each time
tC<-rnorm(1000,mean=15,sd=2)
hist(tC, xlab="size",main="Time C")
text(x = -1, y = 240, labels = "(C)", xpd = NA)
abline(v = mean(tC), col="red", lwd=3, lty=2)

set.seed(89) #to ensure repeatability; skip to vary results slightly each time
tD<-rnorm(1000,mean=16,sd=2)
hist(tD, xlab="size",main="Time D")
text(x = 3, y = 240, labels = "(D)", xpd = NA)
abline(v = mean(tD), col="red", lwd=3, lty=2)

plot(tB~tA, xlab = "size at Time A", ylab = "size at Time B")
text(x = -0.5, y = 22, labels = "(E)", xpd = NA)
summary(lm(tB~tA)) #UNCORRELATED
plot(tC~tB, xlab = "size at Time B", ylab = "size at Time C")
text(x = 3, y = 24, labels = "(F)", xpd = NA)
summary(lm(tC~tB)) #UNCORRELATED
plot(tD~tC, xlab = "size at Time C", ylab = "size at Time D")
text(x = 2, y = 25, labels = "(G)", xpd = NA)
summary(lm(tC~tD)) #UNCORRELATED
plot(tD~tA, xlab = "size at Time A", ylab = "size at Time D")
text(x = -0.5, y = 25, labels = "(H)", xpd = NA)
summary(lm(tD~tA)) #UNCORRELATED

grAdec<-tB-tA
hist(grAdec)
mean(grAdec) #nearly 3

grBdec<-tC-tB
hist(grBdec)
mean(grBdec) #nearly 2

grCdec<-tD-tC
hist(grCdec)
mean(grCdec) #nearly 1

grDdec<-tD-tA
hist(grDdec)
mean(grDdec) #for a total of 6


##now combine into dataframe
A<-as.data.frame(cbind("size"=tA,"growthRate"=grAdec,"InitialSize"=tA))
A$time<-rep("A",1000)
A$diff<-rep("A to A",1000)
head(A)

B<-as.data.frame(cbind("size"=tB,"growthRate"=grAdec,"InitialSize"=tA))
B$time<-rep("B",1000)
B$diff<-rep("A to B",1000)
head(B)

C<-as.data.frame(cbind("size"=tC,"growthRate"=grBdec,"InitialSize"=tB))
C$time<-rep("C",1000)
C$diff<-rep("B to C",1000)
head(C)

D<-as.data.frame(cbind("size"=tD,"growthRate"=grCdec,"InitialSize"=tC))
D$time<-rep("D",1000)
D$diff<-rep("C to D",1000)
head(D)


simDec<-as.data.frame(rbind(A,B,C,D))
simDec$stdGR<-(simDec$growthRate/simDec$InitialSize) #calculate standardized growth rate
simDec$time<-factor(simDec$time,ordered = T,levels=c("A","B","C","D")) #making time ordinal

#and generate plots and statistics

#For size over time
DecSize<-ggplot(simDec, aes(x=time, y=size)) + geom_boxplot() +theme_bw()+ylab(NULL)+xlab(NULL) +ggtitle("Decr. Growth Rate")+ geom_signif(comparisons = list(c("A", "B"),c("B", "C"),c("C", "D")), map_signif_level=TRUE) #simulated size data gets bigger over time steps on average, but no correlation between individuals across time-points
Mds<-aov(size~time,simDec)
summary(Mds)
TukeyHSD(Mds,which="time") #size at each time-step sig diff from all others

#For growth over time
subDec1=subset(simDec,time!="A") #3 pairwise differences for 4 time-points
DecGrowth<-ggplot(subDec1, aes(x=diff, y=growthRate)) + geom_boxplot() +theme_bw() +ylab(NULL) +xlab(NULL) + geom_signif(comparisons = list(c("A to B", "B to C"),c("B to C", "C to D")), map_signif_level=TRUE)
Mdgr<-aov(growthRate~diff,subDec1) 
summary(Mdgr) #sig diff over time
TukeyHSD(Mdgr,which="diff") #size at each time-step sig diff from all others, decreases


#For standardized growth over time
DecGrowthStd<-ggplot(subDec1, aes(x=diff, y=stdGR)) + geom_boxplot() +theme_bw() +ylab(NULL) +xlab(NULL) + geom_signif(comparisons = list(c("A to B", "B to C"),c("B to C", "C to D")), map_signif_level=TRUE)
Mdgrs<-aov(stdGR~diff,subDec1) 
summary(Mdgrs) #decreasing over time
TukeyHSD(Mdgrs,which="diff") #sig difference bewteen first and last time windows


#for growth rate vs initial size 
DecGrowthVSsize<-ggplot(subDec1, aes(x=InitialSize, y=growthRate)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+theme_bw() +ylab(NULL) +xlab(NULL) +ggtitle("Decreasing Growth Rate")
summary(lm(growthRate~InitialSize,subDec1)) #significant negative relationship


###now, plot standardized growth rates and see how relationship looks

DecStdGrowthVSsize<-ggplot(subDec1, aes(x=InitialSize, y=stdGR)) + geom_point(alpha=0.3) +geom_smooth()+theme_bw()+ylab(NULL) +xlab(NULL)
 #very strong negative (and nonlinear) relationship
par(mfrow=c(2,2))
plot(lm(stdGR~InitialSize,subDec1)) #yep, violates assumptions 


######################################################################################
#now for real world data - from Million et al. 2022
#https://github.com/wyattmillion/Acer_Morphological_Plasticity/blob/main/AcerMorphologyData_Imputed2.csv

dat<-read.csv(file="AcerMorphologyData_Imputed2.csv",stringsAsFactors=TRUE,header=TRUE)
summary(dat) #lots of variables, will restrict to just TLE and reformat as simulated datasets above
names(dat) #want columns: 8,15,22,29,36

datTLE<-dat[,c(8,15,22,29,36)]
summary(datTLE)

dat$growth<-dat$TLE12..cm. - dat$TLE1..cm. #calculate growth
dat$stdGR<-(dat$growth/dat$TLE1..cm.) #calculate standardized growth rate

summary(dat) #367 NAs or missing datapoints total - remove

dat2<-dat[complete.cases(dat$growth), ]


par(mfrow=c(2,5))
hist(dat$T0_TLE, xlab="size",main="Time 0")
text(x = -3, y = 90, labels = "(A)", xpd = NA)
abline(v = mean(tA), col="red", lwd=3, lty=2)

hist(dat$T3_TLE, xlab="size",main="Time 3 months")
text(x = -5, y = 90, labels = "(B)", xpd = NA)
abline(v = mean(tB), col="red", lwd=3, lty=2)

hist(dat$T6_TLE, xlab="size",main="Time 6 months") #becomming more non-normal as we go...breakage
text(x = -8, y = 110, labels = "(C)", xpd = NA)
abline(v = mean(tC), col="red", lwd=3, lty=2)

hist(dat$T9_TLE, xlab="size",main="Time 9 months")
text(x = -8, y = 100, labels = "(D)", xpd = NA)
abline(v = mean(tD), col="red", lwd=3, lty=2)

hist(dat$T12_TLE, xlab="size",main="Time 12 months")
text(x = -8, y = 110, labels = "(E)", xpd = NA)
abline(v = mean(tD), col="red", lwd=3, lty=2)

plot(dat$T3_TLE~dat$T0_TLE, xlab = "size at Time 0", ylab = "size at Time 3 mos")
text(x = -3, y = 40, labels = "(F)", xpd = NA)
summary(lm(dat$T3_TLE~dat$T0_TLE)) #positively correlated
plot(dat$T6_TLE~dat$T3_TLE, xlab = "size at Time 3 mos", ylab = "size at Time 6 mos")
text(x = -5, y = 63, labels = "(G)", xpd = NA)
summary(lm(dat$T6_TLE~dat$T3_TLE)) #positively correlated
plot(dat$T9_TLE~dat$T6_TLE, xlab = "size at Time 6 mos", ylab = "size at Time 9 mos")
text(x = -7, y = 82, labels = "(H)", xpd = NA)
summary(lm(dat$T9_TLE~dat$T6_TLE)) #positively correlated
plot(dat$T12_TLE~dat$T9_TLE, xlab = "size at Time 9 mos", ylab = "size at Time 12 mos")
text(x = -8, y = 180, labels = "(I)", xpd = NA)
summary(lm(dat$T12_TLE~dat$T9_TLE)) #strong positive corr

plot(dat$T12_TLE~dat$T0_TLE, xlab = "size at Time 0", ylab = "size at Time 12 mos")
text(x = -3, y = 180, labels = "(J)", xpd = NA)
summary(lm(dat$T12_TLE~dat$T0_TLE)) #positive, but not as strong

#now calculate growth rates between time steps

grAmil<-dat$T3_TLE-dat$T0_TLE
hist(grAmil)
summary(grAmil) #13 NAs due to missingness/mortality
grAmilCC<-grAmil[complete.cases(grAmil)]
mean(grAmilCC) #1.95
median(grAmilCC) #2.37

grBmil<-dat$T6_TLE-dat$T3_TLE
hist(grBmil)
summary(grBmil) #36 NAs due to missingness/mortality
grBmilCC<-grBmil[complete.cases(grBmil)]
mean(grBmilCC) #2.66
median(grBmilCC) #3.32

grCmil<-dat$T9_TLE-dat$T6_TLE
hist(grCmil)
summary(grCmil) #64 NAs due to missingness/mortality
grCmilCC<-grCmil[complete.cases(grCmil)]
mean(grCmilCC) #4.83
median(grCmilCC) #4.56

grDmil<-dat$T12_TLE-dat$T9_TLE
hist(grDmil)
summary(grDmil) #70 NAs due to missingness/mortality
grDmilCC<-grDmil[complete.cases(grDmil)]
mean(grDmilCC) #11.7
median(grDmilCC) #9.38

grEmil<-dat$T12_TLE-dat$T0_TLE
hist(grEmil)
summary(grEmil) #65 NAs due to missingness/mortality
grEmilCC<-grEmil[complete.cases(grEmil)]
mean(grEmilCC) #22.2
median(grEmilCC) #18.51



##now combine into dataframe
A<-as.data.frame(cbind("size"=dat$T0_TLE,"growthRate"=grAmil,"InitialSize"=dat$T0_TLE,"breaks"=dat$CulumativeBreaks))
A$time<-rep("T0",nrow(A))
A$diff<-rep("T0 to T0",nrow(A))
head(A)
summary(A) #need to kill off NAs
Acc<-A[complete.cases(A$growthRate),]
summary(Acc)

B<-as.data.frame(cbind("size"=dat$T3_TLE,"growthRate"=grAmil,"InitialSize"=dat$T0_TLE,"breaks"=dat$T3_Break))
B$time<-rep("T3",nrow(B))
B$diff<-rep("T0 to T3",nrow(B))
head(B)
summary(B) #need to kill off NAs
Bcc<-B[complete.cases(B$growthRate),]
summary(Bcc)


C<-as.data.frame(cbind("size"=dat$T6_TLE,"growthRate"=grBmil,"InitialSize"=dat$T3_TLE,"breaks"=dat$T6_Break))
C$time<-rep("T6",nrow(C))
C$diff<-rep("T3 to T6",nrow(C))
head(C)
summary(C) #need to kill off NAs
Ccc<-C[complete.cases(C$growthRate),]
summary(Ccc)


D<-as.data.frame(cbind("size"=dat$T9_TLE,"growthRate"=grCmil,"InitialSize"=dat$T6_TLE,"breaks"=dat$T9_Break))
D$time<-rep("T9",nrow(D))
D$diff<-rep("T6 to T9",nrow(D))
head(D)
summary(D) #need to kill off NAs
Dcc<-D[complete.cases(D$growthRate),]
summary(Dcc)


E<-as.data.frame(cbind("size"=dat$T12_TLE,"growthRate"=grDmil,"InitialSize"=dat$T9_TLE,"breaks"=dat$T12_Break))
E$time<-rep("T12",nrow(E))
E$diff<-rep("T9 to T12",nrow(E))
head(E)
summary(E) #need to kill off NAs
Ecc<-E[complete.cases(E$growthRate),]
summary(Ecc)



simMil<-as.data.frame(rbind(Acc,Bcc,Ccc,Dcc,Ecc))

summary(simMil) 
head(simMil)

simMil$stdGR<-(simMil$growthRate/simMil$InitialSize) #calculate standardized growth rate - errors arise bc in some instances initial size is 0 due to breakage...must also remove these NAs and infinite values

simMil[!is.finite(simMil$stdGR),] <- NA

summary(simMil) #convert INF to NA

simMil<-simMil[complete.cases(simMil$stdGR),] #now remove NA

simMil$time<-factor(simMil$time,ordered = T,levels=c("T0","T3","T6","T9","T12")) #making time ordinal

subMilnb=subset(simMil,breaks<1) #remove breakage


#and generate plots and statistics

#For size over time
MilSize<-ggplot(subMilnb, aes(x=time, y=size)) + geom_boxplot() +theme_bw()+ylab(NULL)+xlab(NULL) +ggtitle("A. cervicornis")+ geom_signif(comparisons = list(c("T0", "T3"),c("T9", "T12")), map_signif_level=TRUE) #size increases over time 
#NOTE tukey stars are overestimated - need to correct for stats below

Mms<-aov(size~time,subMilnb)
summary(Mms)
TukeyHSD(Mms,which="time") #T0 to T3 not diff, but all other time steps increase

#For growth over time
subMil1=subset(subMilnb,time!="T0") #4 pairwise differences for 5 time-points
MilGrowth<-ggplot(subMil1, aes(x=diff, y=growthRate)) + geom_boxplot() +theme_bw() +ylab(NULL) +xlab(NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_signif(comparisons = list(c("T6 to T9", "T9 to T12")), map_signif_level=TRUE)
Mmgr<-aov(growthRate~diff,subMil1) 
summary(Mmgr) #sig diff over time
TukeyHSD(Mmgr,which="diff") #size at each time-step sig diff from all others, decreases


#For standardized growth over time
MilGrowthStd<-ggplot(subMil1, aes(x=diff, y=stdGR)) + geom_boxplot() +theme_bw() +ylab(NULL) +xlab(NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_signif(comparisons = list(c("T6 to T9", "T9 to T12")), map_signif_level=TRUE)
Mmgrs<-aov(stdGR~diff,subMil1) 
summary(Mmgrs) #decreasing over time
TukeyHSD(Mmgrs,which="diff") #sig difference bewteen first and last time windows


#for growth rate vs initial size, remove instances of breakage
MilGrowthVSsize<-ggplot(subMil1, aes(x=InitialSize, y=growthRate)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+theme_bw() +ylab(NULL) +xlab(NULL) +ggtitle("A. cervicornis")
summary(lm(growthRate~InitialSize,subMil1)) #significant positive relationship
par(mfrow=c(2,2))
plot(lm(stdGR~InitialSize,subMil1)) #violates assumptions bc variance increases over time



###now, plot standardized growth rates and see how relationship looks

MilStdGrowthVSsize<-ggplot(subMil1, aes(x=InitialSize, y=stdGR)) + geom_point(alpha=0.3) +geom_smooth()+theme_bw()+ylab(NULL) +xlab(NULL)
 #very strong negative (and nonlinear) relationship
par(mfrow=c(2,2))
plot(lm(stdGR~InitialSize,subDec1)) #yep, violates assumptions 

head(dat)



##########################################
## one plot to rule them all
###########################################

#Figure 1
pdf(file="sizeGrowthPanelsFig1.pdf", width=8, height=8)
ggarrange(StatSize,IncSize,DecSize,MilSize,StatGrowth,IncGrowth,DecGrowth,MilGrowth,StatGrowthStd,IncGrowthStd,DecGrowthStd,MilGrowthStd,ncol=4,nrow=3,labels=c('A','B','C','D','E','F','G','H','I', 'J', 'K','L'),font.label = list(size = 12, color = "black", face = "plain", family = NULL))
dev.off()

#######
#REGRESSION TO THE MEAN. the negative relationship is a statistical artifact of initial size being included in both X and Y terms 
#regardless of underlying relationship


### NOW lets apply correction for regression to mean using function from Kelley & Price 2005, 
## as scripted by Alex Gunderson for plasticity and modified here for calculating growth
## modifications are simply to change references to plasticity to growth

#R script for function to correct regression values for regression to the mean from Kelly and Price (2005). 
#Takes two vectors, one for each size measurement used to calculate growth, 
#and returns a data frame with three columns: 
#raw.growth – unadjusted growth values. 
#dstar – the correction factors to be applied to the mean raw growth value to create adjusted growth values. 
#adj.growth – adjusted growth values.

#for our purposes, m1 will be size at time A and m2 will be size at time D


#first, read in function - execute all code between ##### below
######################
rttm.adj<-function(m1, m2){
  raw.growth<-m2-m1
  vart<-var.test(m1,m2,paired = T) ## variances equal? 
  vpv<-vart$p.value # var.test p value
  m1m2cor<-cor.test(m1, m2) # test correlation between m1 and m2 
  rho<-m1m2cor$estimate # correlation coefficient between m1 and m2 
  m1sd<-sd(m1) # m1 sd
  m2sd<-sd(m2) # m2 sd
  m1v<-var(m1) # m1 var
  m2v<-var(m2) # m2 var
  m1m<-mean(m1) # m1 mean
  m2m<-mean(m2) # m2 mean
  pm<-mean(raw.growth)
  rho2<-(2*rho*m1sd*m2sd)/(m1v+m2v) # adjusted correlation coefficient used if variances are equal
  rhof<-ifelse(vpv <= 0.05, rho, rho2) # which rho is used for dstar calculation is based on variance comparison
  dstar<-(rhof*(m1-m1m)-(m2-m2m))*-1 # adjustment values. Multiply by -1 to flip sign because Kelly and Price based on plasticity as m1-m2, not m2-m1 as in most thermal tolerance estimates
  adj.growth <- pm+dstar # corrected plasticity. 
  out<-as.data.frame(cbind(raw.growth, dstar, adj.growth)) 
  return(out)
}

##########################

#our focal datasets are:
#subStat1
#subInc1
#subDec1

#where m1 or size at time A = InitialSize
#and m2 or size at time D = size

########### Static growth rate

statAdj<-rttm.adj(subStat1$InitialSize,subStat1$size)

subStat1$growthRate.rttm.adj<-statAdj$adj.growth

#statAdjstdGR<-rttm.adj(subStat1$InitialSize,subStat1$stdGR)

#subStat1$stdGR.rttm.adj<-statAdjstdGR$adj.growth

StatGrowthVSsizeADJ<-ggplot(subStat1, aes(x=InitialSize, y=growthRate.rttm.adj)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+stat_cor(method='pearson')+theme_bw() +ylab("Adjusted Growth") +xlab("Initial Size")
cor(subStat1$InitialSize,subStat1$growthRate.rttm.adj, method='pearson') #no relationship ####### CORRECT BELOW


#StatStdGrowthVSsizeADJ<-ggplot(subStat1, aes(x=InitialSize, y=stdGR.rttm.adj)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+theme_bw() +ylab("Std Growth Adjusted for Regression to Mean") +xlab("Initial Size")
#summary(lm(stdGR.rttm.adj~InitialSize,subStat1)) #spurious positive relationship - cannot apply this correction function in this case!


########### Increasing growth rate
subInc1
incAdj<-rttm.adj(subInc1$InitialSize,subInc1$size)

subInc1$growthRate.rttm.adj<-incAdj$adj.growth

incAdjstdGR<-rttm.adj(subInc1$InitialSize,subInc1$stdGR)

subInc1$stdGR.rttm.adj<-incAdjstdGR$adj.growth


IncGrowthVSsizeADJ<-ggplot(subInc1, aes(x=InitialSize, y=growthRate.rttm.adj)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+stat_cor(method='pearson')+theme_bw() +ylab(NULL) +xlab("Initial Size")
cor(subInc1$InitialSize, subInc1$growthRate.rttm.adj,method='pearson') #positive correlation recovered! although not strong
cor(subInc1$InitialSize, subInc1$growthRate)
########### Decreasing growth rate

decAdj<-rttm.adj(subDec1$InitialSize,subDec1$size)

subDec1$growthRate.rttm.adj<-decAdj$adj.growth

DecGrowthVSsizeADJ<-ggplot(subDec1, aes(x=InitialSize, y=growthRate.rttm.adj)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+stat_cor(method='pearson')+theme_bw() +ylab(NULL) +xlab("Initial Size")
cor(subDec1$InitialSize, subDec1$growthRate.rttm.adj,method='pearson') #negative relationship remains!


########### Million 2022

milAdj<-rttm.adj(subMil1$InitialSize,subMil1$size)

subMil1$growthRate.rttm.adj<-milAdj$adj.growth

MilGrowthVSsizeADJ<-ggplot(subMil1, aes(x=InitialSize, y=growthRate.rttm.adj)) + geom_point(alpha=0.3) +geom_smooth(method='lm')+stat_cor(method='pearson')+theme_bw() +ylab(NULL) +xlab("Initial Size")
cor(subMil1$InitialSize, subMil1$growthRate.rttm.adj,method='pearson') #positive correlation remains - this time strong


#Figure 2
ggarrange(StatGrowthVSsize,IncGrowthVSsize,DecGrowthVSsize,MilGrowthVSsize,StatStdGrowthVSsize,IncStdGrowthVSsize,DecStdGrowthVSsize,MilStdGrowthVSsize,StatGrowthVSsizeADJ,IncGrowthVSsizeADJ,DecGrowthVSsizeADJ,MilGrowthVSsizeADJ,ncol=4,nrow=3,labels=c('A','B','C','D','E','F','G','H','I','J','K','L'),font.label = list(size = 12, color = "black", face = "plain", family = NULL))


##### Now correct for slope for Million 2022 using repeat T0 measures data from Million et al 2021

setwd("~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/RegressionToMean")

#data downloaded from supplement: https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2021.646475/full
#note that asterisks in raw sheet have been replaced with NA as assumption is that this is missing data, then export as CSV

dat2<-read.csv(file="Acer3DMorphologyData.csv",stringsAsFactors=TRUE,header=TRUE)
summary(dat2)

#Blomvquist's formula uses variance in measurement error to calculate true slope

#SlopeCorrected = (SlopeRaw + k) / (1-k) 

#where k = variance of measurement error / total variance

#columns we care about are T0_TLE and T0_FieldTLE (repeat independent measures of TLE to calculate error)

#First calculate measurement error
dat2$meanDif<-dat2$T0_TLE-dat2$T0_FieldTLE

#then calculate variance of measurement error and store as k numerator
k_num<-var(dat2$meanDif)

#calculate total variance in TLE from original dataset
k_den<-var(subMilnb$size)

#then calculate full k
k<-k_num/k_den 
k #pretty small value (expected, as shown in Million 2021 the 2 measures are pretty close)

#now get slope 
lm1<-lm(growthRate.rttm.adj~InitialSize,subMil1)
summary(lm1) #significant positive relationship

SlopeRaw<-coef(lm1)[2] #pull the slope

(SlopeRaw + k) / (1 - k) #correct for measurement error - again, slight adjustment, still positive!

