

#' Hydrolic function for prediction of soil water content
#'
#' @param Wt Soil water content at time t
#' @param Tt Prediction of temperature between time t and t+1
#' @param Pt Prediction of precipitation between time t and t+1
#'
#' @return Predition of soil water content at time t+1
#' @export
Hydro<-function(Wt,Tt,Pt){
  ET= param$hydroETa + param$hydroETb * param$hydroETx*(0.46*Tt + 8.13);
  #if(Pt<param$dryDayTh) f=0 else f= param$hydroKs*( 1-(param$hydroFi*(param$hydroWatS-Wt)/Pt) )
  f=Pt*(1- pow((Wt-param$hydroWatR)/(param$hydroWatS-param$hydroWatR),param$hydroM) );
  g= param$hydroKs*pow((Wt-param$hydroWatR)/(param$hydroWatS-param$hydroWatR),3+2/param$hydroLamba);
  e= ET*(Wt-param$hydroWatR)/(param$hydroWatS-param$hydroWatR);
  return(Wt+f-e-g);
}


#' Set the parameters of the Gaussian SSM.
#'
#' @param t Time period of observing sensor data
#' @param FF design matrix of system equation.
#' @param GG design matrix of observation equation.
#' @param V Observation variance of the nGSSM
#' @param W System variance of the nGSSM (optional values)
#' @param m0 Initial mean of posterior distribution at the insertion time (t=0)
#' @param C0 Initial variance of posterior distribution at the insertion time (t=0)
#'
#' @return A list containing the parameters of the GSSM.
#'
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
#' @export
setModel<-function(param,t,GG,V,W,m0,C0){
  model<-list(t=param$tMax)
  #model$FF<-matrix(data=c(1,0,0.044,1.549),ncol=2)
  model$GG<-1
  model$W<-param$gSSMW
  model$V<-param$gSSMV
  model$m0<-param$gSSMm0
  model$C0<-param$gSSMc0
  return(model)
}


#' Gaussian SSM filtering
#'
#' @param mod set of parameters needed in the GSSM model.
#' @param D Set of the weight and feed intake data generated from the simulation.
#' @param Tem Temperature data.
#' @param Pre Precipitation data.
#' @param W System variance of the GSSM.
#' @param V Observation variance of the GSSM.
#'
#' @return Updated information of filtering.
#'
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
#' @export
DLMfilter<-function(mod,D,Tem,Pre,W,V){

  # The means of posterior
  L1<-c()
  # The variance matrices of posterior
  L2<-c()
  # The variance matrices of prior
  Rt<-c()
  # Mean values of posterior
  meanPos<-c()
  # Variance values of posterior
  varPos<-c()

  for(i in 1:mod$t){
    #Prior
    if(i==1){
      at<-mod$GG * mod$m0
      meanPos[i]<-at
      Rt[i]<-mod$GG * mod$C0 * mod$GG + W
      FF<-Hydro(D[1],Tem[1],Pre[1])
    }
    else{
      at<-mod$GG * L1[i-1]
      meanPos[i]<-at
      Rt[i]<-mod$GG * L2[i-1] * mod$GG + W
      FF<-Hydro(D[i-1],Tem[i-1],Pre[i-1])
    }
    # One step forcast
    ft<-FF * at
    Qt<- FF * Rt[i] * FF + V

    #Posterior (we see y_t here)
    At<-Rt[i] * FF / Qt
    et<-D[i]-ft
    ct<-Rt[i] - At * Qt * At
    varPos[i]<-ct
    L1[i]<-at + At * et
    L2[i]<-Rt[i] - At * Qt * At

  }
  dlm<-list()
  dlm$L1<-L1
  dlm$L2<-L2
  dlm$Rt<-Rt
  dlm$meanPos<-meanPos
  dlm$varPos<-varPos
  return(dlm)
}


#' Gaussian SSM smoothing
#'
#' @param mod set of parameters needed in the GSSM model.
#' @param fdlm1 output information of filtering the GSSM using the function DLMfilter.
#' @param W System variance of the GSSM.
#'
#' @return Updated information of smoothing
#'
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
#' @export
Smoother<-function(mod,fdlm1,W){
  mts<-array(NA,dim=mod$t)
  Cts<-array(NA,dim=mod$t)
  Bt<-array(NA,dim=mod$t)

  mt<-fdlm1$L1
  Ct<-fdlm1$L2
  Rt<-fdlm1$Rt

  mts[mod$t]<-mt[mod$t]
  Cts[mod$t]<-Ct[mod$t]

  for(i in ((mod$t-1):1)){
    Bt[i]<-Ct[i] * mod$GG / Rt[i+1]
    mts[i]<-mt[i] + Bt[i] * mts[i+1] - mod$GG * mt[i]
    Cts[i]<-Ct[i] + Bt[i] * (Cts[i+1] - Rt[i+1]) * Bt[i]
  }

  #for t=0
  Bt0<-mod$C0 * mod$GG / (mod$GG * mod$C0 * mod$GG + W)
  mts0<-mod$m0 + Bt0 * (mts[1] - mod$GG * mod$m0)
  Cts0<-mod$C0 + Bt0 * (Cts[1] - Rt[1]) * t(Bt0)

  smo<-list()
  smo$mts<-mts
  smo$Cts<-Cts
  smo$mts0<-mts0
  smo$Cts0<-Cts0

  return(smo)
}



#' EM algorithm for estimation of system variance in Gaussian SSM
#'
#' @param mod Set of parameters needed in the GSSM model.
#' @param fdlm1 Output information of filtering the GSSM using the function DLMfilter.
#' @param sdlm1 Output information of smoothing the GSSM using the function Smoother.
#' @param D Set of the weight and feed intake data generated from the simulation.
#' @param Wm Updated system variance of the GSSM.
#'
#' @return Updated information of system variance W
#'
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
#' @export
EM<-function(mod,fdlm1,sdlm1,D,Wm){
  Bt<-array(NA,dim=mod$t)
  Lt<-array(NA,dim=mod$t)
  W1<-0
  V1<-0
  result<-list()

  mt<-fdlm1$L1
  Ct<-fdlm1$L2
  mts<-sdlm1$mts
  Cts<-sdlm1$Cts
  mts0<-sdlm1$mts0
  Cts0<-sdlm1$Cts0


  for(i in 1:mod$t){
    if(i==1){
      Bt0<-mod$C0 * mod$GG / (mod$GG * mod$C0 * mod$GG + Wm)
      Lt[i]<-Cts[i] + mod$GG * Cts0 * mod$GG - Cts[i] * Bt0 - Bt0 * Cts[i]
      W1<-Lt[i]+(mts[i] - mod$GG * mts0) * (mts[i] - mod$GG * mts0) + W1
      # FF<-Hydro(D[1],Tem[1],Pre[1])
      # V1<-FF * Cts[i] * FF + (D[i] - FF * mts[i]) * (D[i] - FF * mts[i]) + V1
    }else{
      Bt[i-1]<-Ct[i-1] * mod$GG / (mod$GG * Ct[i-1] * mod$GG + Wm)
      Lt[i]<-Cts[i] + mod$GG * Cts[i-1] * mod$GG - Cts[i] * Bt[i-1] - Bt[i-1] * Cts[i]
      W1<-Lt[i]+(mts[i] - mod$GG * mts[i-1]) * (mts[i] - mod$GG * mts[i-1]) + W1
      # FF<-Hydro(D[i-1],Tem[i],Pre[i])
      # V1<-FF * Cts[i] * FF + (D[i] - FF * mts[i]) * (D[i] - FF * mts[i]) + V1
    }
  }
  W1<-(W1+W1)/2
  # V1<-(V1+V1)/2
  W1<-W1/mod$t
  # V1<-V1/mod$t
  result[[1]]<-W1
  # result[[2]]<-V1
  return(result)
}





#' Simulate the observed and true soil water content given weather data and initial soil-water condition (at time t=1)
#'
#' @param iniTrueWat Initial true soil water content
#' @param temData Weather data regarding temperature
#' @param precData Weather data regarding precipitation
#' @param param Parameter values given in R function \code{setParam}
#'
#' @return A list containing simulated observed and true soil water content.
#' @export
simWat<-function(temData,precData,iniTrueWat,param){
  soilWatTrue<-c()
  soilWatObs<-c()
  soilWatTrue[1]<-iniTrueWat
  soilWatObs[1]<-iniTrueWat

  for(i in 2:param$tMax){
    if(soilWatTrue[i-1]>param$hydroWatS) soilWatTrue[i-1] = param$hydroWatS-1;
    soilWatTrue[i]<- Hydro(soilWatTrue[i-1],temData[i-1],precData[i-1])
    soilWatObs[i]<- soilWatTrue[i]  + mean(rnorm(n = 10, mean=0, sd=sqrt(param$gSSMV)))
  }

  watInfo<-list()
  watInfo$soilWatTrue<-soilWatTrue
  watInfo$soilWatObs<-soilWatObs
  return(watInfo)
}



#' Simulate the observed and true soil water content given weather data and initial soil-water condition (at time t=1)
#'
#' @param param Parameter values given in R function \code{setParam}
#' @param coefPre Coeficiant regarding the change in precipitation
#' @param coefTem Coeficiant regarding the change in temperature
#' @param rndValues A set of random values between 0 and 1 for specifying a wet and rainy day.
#'
#' @return A list containing simulated observed and true soil water content.
#' @export
simWeather<-function(param,coefPre,coefTem,rndValues){
  temData<-c()
  precData<-c()

  #probWet= ( (param$prDryWet)/(1+param$prDryWet-param$prWetWet) )*(1+coefPre)
  probWet= coefPre
  #rndValues<-runif(n = param$tMax, min = 0, max = 1)

  for (i in 1:param$tMax){
    if ( probWet > rndValues[i]  ){
      precData[i] = round(rgamma(n = 1,shape = param$precShape, scale = param$precScale),1)
      while(precData[i]<2) precData[i] = round(rgamma(n = 1,shape = param$precShape, scale = param$precScale),1)
      while(precData[i]>18) precData[i] = round(rgamma(n = 1,shape = param$precShape, scale = param$precScale),1)
    }
    else{
      precData[i]=0
    }


    if (precData[i] <0.25){
      temData[i] = round(rnorm(n = 1, mean = coefTem, sd = sqrt(2)),2)
    }
    else{
      temData[i] = round(rnorm(n = 1, mean = coefTem, sd = sqrt(2)),2)
    }

    # if (precData[i] <0.25){
    #   temData[i] = round(rnorm(n = 1, mean = param$temMeanDry*(1-coefTem), sd = sqrt(2)),2)
    # }
    # else{
    #   temData[i] = round(rnorm(n = 1, mean = param$temMeanWet*(1-coefTem), sd = sqrt(2)),2)
    # }

  }

  wedInfo<-list()
  wedInfo$temData<-temData
  wedInfo$precData<-precData
  return(wedInfo)
}



#' van Genuchten equation for estimating matric water potential given volumetric soil-water conten.
#' Van Genuchten equation can be found in \url{http://www.sciencedirect.com/science/article/pii/S0167198700001549}.
#'
#' @param Wt soil water content in volumetric unit (cm^3/cm^3)
#'
#' @return Matric water potential in hPa unit
#' @export
VanGe<-function(Wt){

  WatR<-param$hydroWatR
  WatS<-param$hydroWatS

  X=pow((Wt-WatR)/(WatS-WatR),1/param$hydroVanM)
  Y=pow((1-X)/X,1/param$hydroVanN )
  H=Y/param$hydroVanAlpha # head pressure in cm
  Z=(H/100)*1000*9.8 # water potential in unit Pa calculated by Z=H*p_w*g (p_w=1000 kg/m^3 and g=9.8 m/s^2)
  return(Z/100)
}


#' Calculate the stree loaded by machine on the field.
#'
#'
#' @export
stressMa<-function(presure, load, depth, matric){

  a<-(1/depth)*pow( (load/pi)*2*(presure/100),0.5)
  b<-atan(a)
  nu<-(7/180000)*pow(matric,2) - (31/1200)*matric + (2590/360)
  if(nu<4) nu=4
  if(nu>6) nu=6
  stress= 2*presure*(1-pow(cos(b),nu) )

}

#' power function
#'
#' @param x A double value
#' @parem y A double value
#'
#' @return a double value for power x^y
#' @export
pow<-function(x,y){
  return(x^y)
}














