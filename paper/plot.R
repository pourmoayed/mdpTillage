
library(reshape2)
library(ggplot2)
library(plyr)
library(grid)
library(tikzDevice)
# # Set MDP parameters for using the basic parameters in the functions
param<-setParam()
# Plot the simulated data of weather info and sensor data with the realts of Gussian SSM for three scenarios
####################################################################################################################################

#Adjust the scaling factor for MP and SP
datS1$MP<-datS1$MP*10
datS2$MP<-datS2$MP*10
datS3$MP<-datS3$MP*10
datS1$SP<-datS1$SP*10
datS2$SP<-datS2$SP*10
datS3$SP<-datS3$SP*10


dtPlot<-rbindlist(list(datS1,datS2,datS3))

dat<-melt(data.frame(dtPlot),
          # ID variables - all the variables to keep but not split apart on
          id.vars=c("scenario","t"),
          # The source columns
          #          measure.vars=c("meanWeight","sdWeight","pricePig", "priceFeed", "pricePiglet", "SP", "SF", "SPi"),
          measure.vars=c("MW","MP","SP","Tem","Pre"),
          # Name of the destination column that will identify the original
          # column that the measurement came from
          variable.name="name",
          value.name="y"
)

dat$scenario<-factor(dat$scenario, labels=c("Scenario 1","Scenario 2","Scenario 3"))
dat$name<-mapvalues(dat$name, from = c("MW","MP","SP","Tem","Pre"),
                    to = c("$ \\hat{\\mu}_t $ \\small($\\frac{cm^3}{cm^3}$)",
                           "$ \\hat{m}_t $ \\small(10)",
                           "$ \\hat{c}_t $ \\small(10)",
                           "$ \\hat{q}_t $ \\small($^{\\circ}{\\rm C}$)",
                           "$ \\hat{p}_t $ \\small(mm)") )

#Plot the data related to the simulation and the SSMs:
tikz("ScenariosPaper_plot.tex", width = 10, height = 7, standAlone=T)

plot<-ggplot(data=dat, aes(x=factor(t), y=y, group=name, shape=name, linetype=name ) ) +
  geom_line() + scale_y_continuous(breaks=seq(0,50,1), labels = c(0:50) ) +
  #geom_point() +
  facet_grid(. ~ scenario) +
  xlab("Day number") + ylab(" ")
g <- guide_legend("",nrow=1,byrow=TRUE, override.aes = list(fill=NA))

plot + guides(shape = g, linetype=g)  +
  #geom_histogram(stat="identity", data=datPigs, alpha = 1/4, colour=NA, width=0.25) +
  #  geom_vline(aes(xintercept = w), data=vline.fm, color="gray") +
  #  geom_vline(aes(xintercept = w), data=vline.th, color="gray", linetype="twodash") +
  geom_line() +
  #  geom_text(data=datMarketing, mapping=aes(x=t, y=-1, label=optimalLable), size=4) +
  #  scale_y_discrete(breaks = -2:15, labels = c("","Optimal decisions",0:15) ) +
  theme_bw() +
  theme(legend.position="bottom", panel.background = element_blank(),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.key.width = unit(2, "cm"), legend.text.align=0.5, axis.title.x= element_text(vjust = -0.7),
        axis.text.x = element_text(size=5), axis.text.y = element_text(size=7),
        strip.background=element_rect(fill = NA))

dev.off()

####################################################################################################################################




# Plot optimal tillage decisions for three scenarios
####################################################################################################################################
# Plot optimal decisions
tikz("OptimalPaper_plot.tex", width = 9, height = 6, standAlone=T)

plot(c(0,31), c(1,12), yaxt="n", xlab="Day number", ylab='', xaxt="n", bty='n', pch=NA)
abline(v=1:30,lty= 2, col="gray85")
axis(1, las=1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30), labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"),
     cex.axis=0.7)
axis(2, las=3, col = "white", at=c(2, 5, 8, 11), labels=c("Ploughing", "Harrowing-1", "Harrowing-2", "Planting"), line=1, cex.axis=0.9)
axis(2, las=1, col = "white", at=c(1.5, 2, 2.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-3, cex.axis=0.7)
axis(2, las=1, col = "white", at=c(4.5, 5, 5.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-3, cex.axis=0.7)
axis(2, las=1, col = "white",  at=c(7.5, 8, 8.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-3, cex.axis=0.7)
axis(2, las=1, col = "white", at=c(10.5, 11, 11.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-3, cex.axis=0.7)


#title(main="Optimal ", cex.main=1)

dat<-datS1
dat<-as.data.table(dat)
vecFCor<-c(1.5, 4.5, 7.5, 10.5)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="black", cex = 1.85)
  }
}


dat<-datS2
dat<-as.data.table(dat)
vecFCor<-c(2, 5, 8, 11)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="green", cex = 1.85)
  }
}



dat<-datS3
dat<-as.data.table(dat)
vecFCor<-c(2.5, 5.5, 8.5, 11.5)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="blue", cex = 1.85)
  }
}

dev.off()
####################################################################################################################################




# plot optimal tilaage decisions for Group 1 and 2 in the paper with different weights of trafficability, workability and completion criteria
#####################################################################################################################################

tikz("ComparePolicyPaper_plot.tex", width = 10, height = 5.5, standAlone=T)

par(mar = c(4, 4, 3, 0))
layout(matrix(c(1,1,2,2), 2, 2, byrow = F))

#read policy for Group 1
policy <- read.csv2("polices/weight_high_traf/policyMDP.csv", stringsAsFactors = F)

coefPre1<-0.2
coefPre2<-0.4
coefPre3<-0.6

coefTem1<-16
coefTem2<-14
coefTem3<-12

set.seed(10000)
rndValues<-runif(n = param$tMax, min = 0, max = 1)
#length(rndValues[rndValues<0.25])

iniTrueWat<-30
givenWatInfo<-FALSE

dataOptimal1<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre1,coefTem1,rndValues)$temData, precData = simWeather(param,coefPre1,coefTem1,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )
dataOptimal2<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre2,coefTem2,rndValues)$temData, precData = simWeather(param,coefPre2,coefTem2,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )
dataOptimal3<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre3,coefTem3,rndValues)$temData, precData = simWeather(param,coefPre3,coefTem3,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )

datS1 <- dataOptimal1
datS2 <- dataOptimal2
datS3 <- dataOptimal3
datS1$scenario<-1
datS2$scenario<-2
datS3$scenario<-3
datS1$MP<-datS1$MP*10
datS2$MP<-datS2$MP*10
datS3$MP<-datS3$MP*10
datS1$SP<-datS1$SP*10
datS2$SP<-datS2$SP*10
datS3$SP<-datS3$SP*10

#Plot Group 1
plot(c(0,31), c(1,12), yaxt="n", xlab="Day number", ylab='', xaxt="n", bty='n', pch=NA)
abline(v=1:30,lty= 2, col="gray85")
axis(1, las=1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30), labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"),
     cex.axis=0.7)
title(main="Group 1: $ \\alpha^{\\mathtt{\\small traf}}=0.8, \\alpha^{\\mathtt{\\small work}}=0.2, \\alpha^{\\mathtt{\\small comp}}=1 $", cex.main=1)
axis(2, las=3, col = "white", at=c(2, 5, 8, 11), labels=c("Ploughing", "Harrowing-1", "Harrowing-2", "Planting"), line=1, cex.axis=0.9)
axis(2, las=1, col = "white", at=c(1.5, 2, 2.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)
axis(2, las=1, col = "white", at=c(4.5, 5, 5.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)
axis(2, las=1, col = "white",  at=c(7.5, 8, 8.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)
axis(2, las=1, col = "white", at=c(10.5, 11, 11.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)

dat<-datS1
dat<-as.data.table(dat)
vecFCor<-c(1.5, 4.5, 7.5, 10.5)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="black", cex = 1.5)
  }
}


dat<-datS2
dat<-as.data.table(dat)
vecFCor<-c(2, 5, 8, 11)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="green", cex = 1.5)
  }
}



dat<-datS3
dat<-as.data.table(dat)
vecFCor<-c(2.5, 5.5, 8.5, 11.5)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="blue", cex = 1.5)
  }
}


#read policy for Group 2
policy <- read.csv2("polices/weight_high_work/policyMDP.csv", stringsAsFactors = F)

coefPre1<-0.2
coefPre2<-0.4
coefPre3<-0.6

coefTem1<-16
coefTem2<-14
coefTem3<-12

set.seed(10000)
rndValues<-runif(n = param$tMax, min = 0, max = 1)
#length(rndValues[rndValues<0.25])

iniTrueWat<-30
givenWatInfo<-FALSE

dataOptimal1<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre1,coefTem1,rndValues)$temData, precData = simWeather(param,coefPre1,coefTem1,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )
dataOptimal2<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre2,coefTem2,rndValues)$temData, precData = simWeather(param,coefPre2,coefTem2,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )
dataOptimal3<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre3,coefTem3,rndValues)$temData, precData = simWeather(param,coefPre3,coefTem3,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )

datS1 <- dataOptimal1
datS2 <- dataOptimal2
datS3 <- dataOptimal3
datS1$scenario<-1
datS2$scenario<-2
datS3$scenario<-3
datS1$MP<-datS1$MP*10
datS2$MP<-datS2$MP*10
datS3$MP<-datS3$MP*10
datS1$SP<-datS1$SP*10
datS2$SP<-datS2$SP*10
datS3$SP<-datS3$SP*10



#plot the Group 2
plot(c(0,31), c(1,12), yaxt="n", xlab="Day number", ylab='', xaxt="n", bty='n', pch=NA)
abline(v=1:30,lty= 2, col="gray85")
axis(1, las=1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30), labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"),
     cex.axis=0.7)
title(main="Group 2: $ \\alpha^{\\mathtt{\\small traf}}=0.2, \\alpha^{\\mathtt{\\small work}}=0.8, \\alpha^{\\mathtt{\\small comp}}=1 $", cex.main=1)
axis(2, las=3, col = "white", at=c(2, 5, 8, 11), labels=c("Ploughing", "Harrowing-1", "Harrowing-2", "Planting"), line=1, cex.axis=0.9)
axis(2, las=1, col = "white", at=c(1.5, 2, 2.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)
axis(2, las=1, col = "white", at=c(4.5, 5, 5.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)
axis(2, las=1, col = "white",  at=c(7.5, 8, 8.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)
axis(2, las=1, col = "white", at=c(10.5, 11, 11.5), labels=c("Scenario 1", "Scenario 2", "Scenario 3"), line=-2.5, cex.axis=0.7)

dat<-datS1
dat<-as.data.table(dat)
vecFCor<-c(1.5, 4.5, 7.5, 10.5)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="black", cex = 1.5)
  }
}


dat<-datS2
dat<-as.data.table(dat)
vecFCor<-c(2, 5, 8, 11)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="green", cex = 1.5)
  }
}



dat<-datS3
dat<-as.data.table(dat)
vecFCor<-c(2.5, 5.5, 8.5, 11.5)
for(j in 1:param$opNum){
  for(i in dat[dat$operation==j,]$t){

    if( (dat[dat$t==i]$optAction=="do.") || (dat[dat$t==i]$optAction=="doF.") ) labCol=15 else labCol=17
    datL<-rep(vecFCor[j],1)
    points(x = c(i), y = datL, type="p", pch=labCol, col="blue", cex = 1.5)
  }
}


dev.off()


#####################################################################################################################################





