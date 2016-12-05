
#' Set the parameters used when build the MDP model
#'
#' @param opNum number of operations for scheduling
#' @param tMax Maximum length of growing cycle.
#' @param opSeq A vector showing operations indexes and sequence
#' @param opE A vector showing earliest start time of operations
#' @param opL A vector showing latest start time of operationss
#' @param opD A vector showing the time needed to complete operation (obtained by machine capacity and the field area)
#' @param opDelay A vector including the delay times after finishing  the operations (except the last operation)
#' @param opFixCost A vector including the fixed costs of tillage operations
#' @param watTh A vector including threshold values used as the optimal level of soil-water content for performing tillage operations (%)
#' @param coefLoss Coeficiant showing yield reduction when soil-water content is not appropriate
#' @param priceYield Price per kg yield (DDK)
#' @param machCap A vector showing the machine capacity for tillage operations
#' @param yieldHa Estimated yield per hetar in field (kg)
#' @param fieldArea Area of field (ha)
#' @param coefTimeliness Timeliness Coeficiant showing yield reduction resulting from postponing a tillage operation.
#' @param costSkip Cost of skipping tillage operations for the current cropping period
#' @param watUpper Upper limits of water content for workability criterion calculated based on the method in \url{http://www.sciencedirect.com/science/article/pii/S0167198700001549}
#' @param watLower Lower limits of water content for workability criterion calculated based on the method in \url{http://www.sciencedirect.com/science/article/pii/S0167198700001549}
#' @param stress Stree of the soil for different values of soil-water content given in \var{centerPointsAvgWat}. The stress values are computed using Terramino (\url{http://www.terranimo.dk/}) and the results are in directory \dir{R/data/stree_strength_machine}
#' @param strength Strength of the soil for different values of soil-water content given in \var{centerPointsAvgWat}. The stress values are computed using Terramino (\url{http://www.terranimo.dk/}) and the results are in directory \dir{R/data/stree_strength_machine}
#' @param weightCompletion Weight of the completion criterion used in the reward function of the MDP.
#' @param weightWorkable Weight of the workability criterion used in the reward function of the MDP.
#' @param weightTraffic Weight of the trafficability criterion used in the reward function of the MDP.
#' @param minOpt The lower tail of the interval related to the best time period for finishing tillage operations.
#' @param maxOpt The upper tail of the interval related to the best time period for finishing tillage operations.
#' @param centerPointsAvgWat Center points for discritization of estimated mean of soil-water content
#' @param centerPointsSdWat Center points for discritization of estimated standard deviation of soil-water content
#' @param centerPointsSdPos Center points for discritization of estimated standard deviation of posterior distribution in SSM
#' @param centerPointsMeanPos Center points for discritization of estimated mean of posterior distribution in SSM
#' @param centerPointsTem Center points for discritization of weather forecast data regarding air temperature
#' @param centerPointsPre Center points for discritization of weather forecast data regarding precipitation
#' @param temMeanDry Mean of air temprature in dry days (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param temMeanWet Mean of air temprature in wet days (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param temVarDry variance of air temprature in dry days (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param temVarWet variance of air temprature in wet days (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param dryDayTh Threshold of precipitation amount for being a dry day (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param precShape Shape parameter of gamma distribution for precipitation amount (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param precScale Scale parameter of gamma distribution for precipitation amount (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param prDryWet Probability of a wet day when the previous day is dry (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param prWetWet Probability of a wet day when the previous day is wet (\url{http://link.springer.com/article/10.1007/BF00142466}).
#' @param hydroWatR Residual soil moisture value (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroWatS Soil moisture value at saturation (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroM Parameter used in the infiltration process of rainfall-runoff model (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroKs Saturated hydralic conductivity used for drainage process in the rainfall-runoff model (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroLamba Parameter used for drainage process in therainfall-runoff model (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroFi Parameter used in the infiltration process of rainfall-runoff model (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroETa Parameter used for evapotranpiration process in the rainfall-runoff model (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroETb Parameter used for evapotranpiration process in the rainfall-runoff model (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroETx Parameter used for evapotranpiration process in the rainfall-runoff model (\url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract}).
#' @param hydroVanAlpha Parameter used in van Genuchten equation to calculate the matric water potensial (\url{http://www.sciencedirect.com/science/article/pii/S0016706198001323})
#' @param hydroVanM Parameter used in van Genuchten equation to calculate the matric water potensial (\url{http://www.sciencedirect.com/science/article/pii/S0016706198001323})
#' @param hydroVanN Parameter used in van Genuchten equation to calculate the matric water potensial (\url{http://www.sciencedirect.com/science/article/pii/S0016706198001323})
#' @param hydroBulkDensity Parameter used to change the unit of soil-water content from kg/kg to cm^3/cm^3(\url{http://www.sciencedirect.com/science/article/pii/S0016706198001323})
#' @param gSSMW System variance of Gaussian SSM.
#' @param gSSMV Observation variance of Gaussian SSM.
#' @param gSSMm0 Initial posterior mean of Gaussian SSM.
#' @param gSSMc0 Initial posterior variance of Gaussian SSM.
#' @param nGSSMm0 Initial posterior mean of non-Gaussian SSM.
#' @param nGSSMc0 Initial posterior mean of non-Gaussian SSM.
#' @param nGSSMK Number of observations in non-Gaussian SSM.
#' @param rewRisk A boolean variable specifing to calculate the reward based on cost parameters or the satisfaction level for trafficability, workability and completion criteria.
#' @param check Check model e.g. do trans pr sum to one
#'
#' @return A list containing all the parameters used in three-level HMDP
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
#' @export
setParam<-function(
  opNum=4,
  tMax=30,
  opSeq=c(1,2,3,4),
  opD=c(6,4,4,4),
  opE=c(1,1+opD[1],1+opD[1]+opD[2],1+opD[1]+opD[2]+opD[3]),#c(1,3,5,7),
  opL=c(tMax-opD[4]-opD[3]-opD[2],tMax-opD[4]-opD[3],tMax-opD[4],tMax),#c(6,8,10,12),
  opDelay=c(0,0,0),
  opFixCost=c(1000,1000,1000,1000),
  watTh=c(50,50,50,50),
  coefLoss=0.3,
  priceYield=10,
  machCap=50,
  yieldHa=100,
  fieldArea=100,
  coefTimeliness=0.01,
  costSkip= 80000,

  watUpper=c(37,37,37,37),
  watLower=c(23.9,23.9,23.9,23.9),
  stress=c(133,133,133,133,117,101,90,78,65),
  strength=c(504,504,504,386,218,138,93,62,35),
  weightCompletion=1,#1,#1,
  weightWorkable=1,#0.8,#0.4,
  weightTraffic=1,#0.2,#0.2,
  minOpt=18,
  maxOpt=23,

  centerPointsAvgWat=seq(hydroWatR+1,hydroWatS, by=5),#seq(5,60, by=6),
  centerPointsSdWat=seq(1,1, by=1),
  centerPointsSdPos=c(0.005,0.01,0.05,0.1),
  centerPointsMeanPos=seq(0.8,1.2,by=0.07),
  centerPointsTem=seq(10,23,by=2.5),
  centerPointsPre=c(0,dryDayTh*2, seq(1,18,by=3) ),

  # centerPointsAvgWat=seq(hydroWatR+1,hydroWatS, by=5),#seq(5,60, by=6),
  # centerPointsSdWat=seq(1,1, by=1),
  # centerPointsSdPos=seq(0.002, 0.009, by=0.03),#c(0.005,0.01,0.05,0.1),#seq(0.005,0.1, by=0.01),
  # centerPointsMeanPos=seq(1,1,by=1),
  # centerPointsTem=seq(10,23,by=3.5),
  # centerPointsPre= c(0,dryDayTh*2, seq(4,15,by=5) ),


  temMeanDry=14.3, #shoule be estimates
  temMeanWet=14.6, #shoule be estimates
  temVarDry=6.2, #shoule be estimates
  temVarWet=4, #shoule be estimates
  dryDayTh=0.25, #shoule be estimates
  precShape=0.36, #shoule be estimates
  precScale=5.1, #shoule be estimates
  prDryWet= 0.176,  #shoule be estimates
  prWetWet=0.71, #shoule be estimates

  # Soil medium
  hydroWatR=0.010*100, #(%)
  hydroWatS=0.439*100, #(%)
  hydroM=15, #shoule be estimates
  hydroKs=12.061, # since it should be converted in mm/day IS IT TRUE?
  hydroLamba=0.085*log(hydroKs)+0.1574, #shoule be estimates
  hydroFi=54.727*log(hydroKs)-323.9, #shoule be estimates
  hydroETa=-2, #shoule be estimates
  hydroETb=1.26, #shoule be estimates
  hydroETx=0.25, #shoule be estimates
  hydroVanAlpha=0.0314,
  hydroVanM=0.1528,
  hydroVanN=1.1804,
  hydroBulkDensity=1.7, # ?

  #Soil fine
  # hydroWatR=0.010*100, #(%)
  # hydroWatS=0.52*100, #(%)
  # hydroM=15, #shoule be estimates
  # hydroKs=24.8, # since it should be converted in mm/day IS IT TRUE?
  # hydroLamba=0.085*log(hydroKs)+0.1574, #shoule be estimates
  # hydroFi=54.727*log(hydroKs)-323.9, #shoule be estimates
  # hydroETa=-2, #shoule be estimates
  # hydroETb=1.26, #shoule be estimates
  # hydroETx=0.25, #shoule be estimates
  # hydroVanAlpha=0.0367,
  # hydroVanM=0.0919,
  # hydroVanN=1.1012,
  # hydroBulkDensity=1.7, # ?


  gSSMW=0.063, #0.01,
  gSSMV=9,
  gSSMm0=1,
  gSSMc0=0.001,
  nGSSMm0=4,
  nGSSMc0=2,
  nGSSMK=20,

  rewRisk=TRUE,

  check = FALSE
){
   model<-list(opNum=opNum)
   model$opSeq<-opSeq
   model$tMax<-tMax
   model$opE<-opE
   model$opL<-opL
   model$opD<-opD
   model$opDelay<-opDelay
   model$opFixCost<-opFixCost
   model$watTh<-watTh
   model$coefLoss<-coefLoss
   model$priceYield<-priceYield
   model$machCap<-machCap
   model$yieldHa<-yieldHa
   model$fieldArea<-fieldArea
   model$coefTimeliness<-coefTimeliness
   model$costSkip<-costSkip

   model$watUpper<-watUpper
   model$watLower<-watLower
   model$stress<-stress
   model$strength<-strength
   model$weightCompletion<-weightCompletion
   model$weightWorkable<-weightWorkable
   model$weightTraffic<-weightTraffic
   model$minOpt<-minOpt
   model$maxOpt<-maxOpt

   model$temMeanDry<-temMeanDry
   model$temMeanWet<-temMeanWet
   model$temVarDry<-temVarDry
   model$temVarWet<-temVarWet
   model$dryDayTh<-dryDayTh
   model$precShape<-precShape
   model$precScale<-precScale
   model$prDryWet<-prDryWet
   model$prWetWet<-prWetWet

   model$hydroWatR<-hydroWatR
   model$hydroWatS<-hydroWatS
   model$hydroM<-hydroM
   model$hydroKs<-hydroKs
   model$hydroFi<-hydroFi
   model$hydroLamba<-hydroLamba
   model$hydroETa<-hydroETa
   model$hydroETb<-hydroETb
   model$hydroETx<-hydroETx

   model$hydroVanAlpha<-hydroVanAlpha
   model$hydroVanM<-hydroVanM
   model$hydroVanN<-hydroVanN
   model$hydroBulkDensity<-hydroBulkDensity

   model$gSSMW<-gSSMW
   model$gSSMV<-gSSMV
   model$gSSMm0<-gSSMm0
   model$gSSMc0<-gSSMc0
   model$nGSSMm0<-nGSSMm0
   model$nGSSMc0<-nGSSMc0
   model$nGSSMK<-nGSSMK

   model$rewRisk <-rewRisk
   model$check <-check

   model$centerPointsAvgWat<-centerPointsAvgWat
   model$centerPointsSdWat<-centerPointsSdWat
   model$centerPointsSdPos<-centerPointsSdPos
   model$centerPointsMeanPos<-centerPointsMeanPos
   model$centerPointsTem<-centerPointsTem
   model$centerPointsPre<-centerPointsPre

   #Discritization of continious states:

   obj<-Discretize()
   disAvgWat<-matrix()
   disSdWat<-matrix()
   disSdPos<-matrix()
   disMeanPos<-matrix()
   disTem<-matrix()
   disPre<-matrix()

   disAvgWat<-as.matrix(obj$discretize1DVec(centerPointsAvgWat, mInf=-1000, inf=1000, asDF=F), ncol=3)
   disSdWat<-as.matrix(obj$discretize1DVec(centerPointsSdWat, inf=100, mInf=0.01, asDF=F), ncol=3)
   disSdPos<-as.matrix(obj$discretize1DVec(centerPointsSdPos, inf=100, mInf=0, asDF=F), ncol=3)
   disMeanPos<-as.matrix(obj$discretize1DVec(centerPointsMeanPos, inf=100, asDF=F), ncol=3)
   disTem<-as.matrix(obj$discretize1DVec(centerPointsTem, inf=100, asDF=F), ncol=3)
   disPre<-as.matrix(obj$discretize1DVec(centerPointsPre, inf=100, mInf=0, asDF=F), ncol=3)

   model$disAvgWat<-disAvgWat
   model$disSdWat<-disSdWat
   model$disSdPos<-disSdPos
   model$disMeanPos<-disMeanPos
   model$disTem<-disTem
   model$disPre<-disPre

   return(model)
}

#' Function to find the index of a state given state value
#'
#' @param st State value
#' @param dis A matrix containing the discritization values of a state variable
#'
#' @return An index related to value st
#' @export
findIndex<-function(st,dis){
  for(i in 1:dim(dis)[1]){
    if( ( st>=dis[i,2] ) & ( st<dis[i,3] )  )
      return(i-1)
  }
  cat("error in index","\n")
  return(-1)
}




