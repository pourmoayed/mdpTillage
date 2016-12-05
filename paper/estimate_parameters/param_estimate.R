
#---------------------------------------------------------------------------------------------------------

# This file contains the codes and functions for reading and analysing the meteorological and sensor data for
# estimation of the parameters for the optimization and statistical models.

# Note: make the current directory ./estimate_parameters as the root directory

#---------------------------------------------------------------------------------------------------------


# Read the soil moisture and weather data for three months in a real farm

metData <- read.table( "sensor_weather_3months.csv", sep = "," , header = T)

metData$X <- NULL


#Calculate average daily data for moisture, temperature, rain, and solar
#----------------------------------------------------------------------------------------------------
days <- seq(from=as.Date('2016-03-4'), to=as.Date("2016-05-1"),by='days' )
metDataDaily <- as.data.frame( matrix(nrow = length(days), ncol = 5)  )
colnames(metDataDaily) <- c("date","meanMois","meanTemp","meanRain", "meanSolar")
for ( i in seq_along(days) ) {
  metDataFilter <- metData[grep(days[i], metData$Timestamp), ]
  meanTemp <- mean(metDataFilter$SLHT5.Air..Air.temperature...C.)
  meanMois <- mean(metDataFilter$EC5..Soil.moisture.....)
  meanRain <- sum(metDataFilter$RainMeter..Rain...mm.)
  meanSolar <- mean(metDataFilter$Solar.panel..Solar.radiation.level...Levels.)
  metDataDaily$date[i] <- paste( days[i] )
  metDataDaily$meanMois[i] <- meanMois
  metDataDaily$meanTemp[i] <- meanTemp
  metDataDaily$meanRain[i] <- meanRain
  metDataDaily$meanSolar[i] <- meanSolar
}
#----------------------------------------------------------------------------------------------------------

#Estimate system variance in the Gaussian SSM using EM algorithm

library(mdpTillage)
require(discretizeGaussian)
# Set HMDP parameters:
param<-setParam()
# Read the data from data base
D<-metDataDaily$meanMois
Tem<-metDataDaily$meanTemp
Pre<-metDataDaily$meanRain
mod<-setModel(param)

# param$hydroKs<-120.61
# param$hydroM<-15
# param$hydroETx<-0.5
# forecast<-c()
# forecast[1]<-30#D[1]
# for(i in 2:length(D)){
#   forecast[i]<-Hydro(forecast[i-1],Tem[i],Pre[i])
# }
# plot(y=forecast,x=1:length(D), type = "l", col = "red", ylim=c(0,70))
# lines(y=D,x=1:length(D),xaxt="n",yaxt="n",xlab="",ylab="", col="blue")


# Implement the EM algorithm for z=1000 iterations :
z<-1000
fdlm<-list()
sdlm<-list()
We<-array(NA,dim=z)
Ve<-array(NA,dim=z)
We1<-0
Ve1<-0

  for(u in 1:z){
    if(u==1){
      fdlm[[u]]<-DLMfilter(mod,D,Tem,Pre,mod$W,mod$V)
      sdlm[[u]]<-Smoother(mod,fdlm[[u]],mod$W)
      We[u]<-EM(mod,fdlm[[u]],sdlm[[u]],D,mod$W)[[1]]
      # Ve[u]<-EM(mod,fdlm[[u]],sdlm[[u]],D,mod$W)[[2]]
    }
    else{
      fdlm[[u]]<-DLMfilter(mod,D,Tem,Pre,We[u-1],mod$V)
      sdlm[[u]]<-Smoother(mod,fdlm[[u]],We[u-1])
      We[u]<-EM(mod,fdlm[[u]],sdlm[[u]],D,We[u-1])[[1]]
      # Ve[u]<-EM(mod,fdlm[[u]],sdlm[[u]],D,We[u-1])[[2]]
    }
  }
  We1<-We[z] + We1
  # Ve1<-Ve[z] + Ve1

#--------------------------------------------------------------------------------------------------------------------------


# Read weather data for september in 2015 and 2016 and estimate parameters for distributions related to weather info.

wedData1 <- read.table( "weatherData_2015_sep", sep = ";" , header = T) # Read data for Sept. 2015
wedData2 <- read.table( "weatherData_2016_sep", sep = ";" , header = T) # Select data for Sept. 2016

# Calculate averge mean and variance of temperature for wed and dry days

avgTem15<-(wedData1$high_temperature + wedData1$low_temperature)/2
avgTem15Pre<-(wedData1[wedData1$precipitation.mm.!=0,]$high_temperature + wedData1[wedData1$precipitation.mm.!=0,]$low_temperature)/2
avgTem15Dry<-(wedData1[wedData1$precipitation.mm.==0,]$high_temperature + wedData1[wedData1$precipitation.mm.==0,]$low_temperature)/2

meanTem15Pre<-mean(avgTem15Pre)
varTem15Pre<-var(avgTem15Pre)
meanTem15Dry<-mean(avgTem15Dry)
varTem15Dry<-var(avgTem15Dry)


avgTem16<-(wedData2$high_temperature + wedData2$low_temperature)/2
avgTem16Pre<-(wedData2[wedData2$precipitation.mm.!=0,]$high_temperature + wedData2[wedData2$precipitation.mm.!=0,]$low_temperature)/2
avgTem16Dry<-(wedData2[wedData2$precipitation.mm.==0,]$high_temperature + wedData2[wedData2$precipitation.mm.==0,]$low_temperature)/2

meanTem16Pre<-mean(avgTem16Pre)
varTem16Pre<-var(avgTem16Pre)
meanTem16Dry<-mean(avgTem16Dry)
varTem16Dry<-var(avgTem16Dry)

meanTemPre<-(meanTem16Pre+meanTem15Pre)/2
varTemPre<-(varTem16Pre+varTem15Pre)/2
meanTemDry<-(meanTem16Dry+meanTem15Dry)/2
varTemDry<-(varTem16Dry+varTem15Dry)/2

# Calculate the probability of wed day following a dry day and probability of wet day following a wet day

counterDry<-0
for(i in 1:(dim(wedData1)[1]) )
  if( (wedData1$precipitation.mm.[i]==0) ) counterDry<-counterDry+1

counterWet<-0
for(i in 1:(dim(wedData1)[1]) )
  if( (wedData1$precipitation.mm.[i]!=0) ) counterWet<-counterWet+1

counterDryWet<-0
for(i in 1:(dim(wedData1)[1]-1) )
  if( (wedData1$precipitation.mm.[i]==0) && (wedData1$precipitation.mm.[i+1]!=0) ) counterDryWet<-counterDryWet+1

counterWetWet<-0
for(i in 1:(dim(wedData1)[1]-1) )
  if( (wedData1$precipitation.mm.[i]!=0) && (wedData1$precipitation.mm.[i+1]!=0) ) counterWetWet<-counterWetWet+1

prDryWet15<-counterDryWet/counterDry
prWetWet15<-counterWetWet/counterWet


counterDry<-0
for(i in 1:(dim(wedData2)[1]) )
  if( (wedData2$precipitation.mm.[i]==0) ) counterDry<-counterDry+1

counterWet<-0
for(i in 1:(dim(wedData2)[1]) )
  if( (wedData2$precipitation.mm.[i]!=0) ) counterWet<-counterWet+1

counterDryWet<-0
for(i in 1:(dim(wedData2)[1]-1) )
  if( (wedData2$precipitation.mm.[i]==0) && (wedData2$precipitation.mm.[i+1]!=0) ) counterDryWet<-counterDryWet+1

counterWetWet<-0
for(i in 1:(dim(wedData2)[1]-1) )
  if( (wedData2$precipitation.mm.[i]!=0) && (wedData2$precipitation.mm.[i+1]!=0) ) counterWetWet<-counterWetWet+1

prDryWet16<-counterDryWet/counterDry
prWetWet16<-counterWetWet/counterWet

prDryWet<-(prDryWet16 +prDryWet15)/2
prWetWet<-(prWetWet16 +prWetWet15)/2


# Calculate the shape and scale of the gamma distribution related to precipitation.

meanPre15<-mean(wedData1$precipitation.mm.)
varPre15<-var(wedData1$precipitation.mm.)
scale15<-varPre15/meanPre15
shape15<-meanPre15/scale15


meanPre16<-mean(wedData2$precipitation.mm.)
varPre16<-var(wedData2$precipitation.mm.)
scale16<-varPre16/meanPre16
shape16<-meanPre16/scale16

meanPre<-(meanPre15+meanPre16)/2
varPre<-(varPre16 + varPre15)/2
scalePre<-varPre/meanPre
shapePre<-meanPre/scalePre

#----------------------------------------------------------------------------------------------------------

#Calculating the upper and lower limit for workability criterion

library(mdpTillage)
# Set MDP parameters:
param<-setParam()

# Based on the paper \url{http://www.sciencedirect.com/science/article/pii/S0167198700001549}
# the optimal limit and upper limit are calculated as follows:

opt<-(param$hydroWatS - param$hydroWatR) * pow(1 + 1/param$hydroVanM, -param$hydroVanM ) + param$hydroWatR
upper<- opt + (param$hydroWatS - opt)*0.4

# In order to calculate the lower limit we use a rule of thumb given in \url{http://www.sciencedirect.com/science/article/pii/S0167198700001549}
# that the water content at lower limit is "the water content at which the strength of the soil is twice the strength at the optimum leve"
# Therefore first using the Van Genuchten equation we calculate the head presure at optimal level of water content and next by using the
# TERRANIMO app. (\url{http://www.terranimo.dk/} ) we calculate the soil strength.

# for soil medium:
VanGe(opt) # head pressure at optimal level (hPa)
# we use this value in TERRANIMO and the soil strength is calculated strengthOpt=89 Kpa and therefore the soil strength  for lowe level
# should be about 180 kPa that using TERRANIMO is equivallent with head presure 1000 hPa . Based on Van Genuchten equation
# the water content 23.9 (cm^3/cm^3)   generates this presure
VanGe(23.9)

# for soil fine:
VanGe(opt) # head pressure at optimal level (hPa)
# we use this value in TERRANIMO and the soil strength is calculated strengthOpt=105 Kpa  and therefore the soil strength  for lowe level
# should be about 210 kPa that using TERRANIMO is equivallent with head presure 1500 hPa . Based on Van Genuchten equation
# the water content 34.9 (cm^3/cm^3)   generates this presure
VanGe(34.9)












