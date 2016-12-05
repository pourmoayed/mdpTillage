


#' Hydrolic function for prediction of soil water content
#'
#' @param param Parameter values given in R function \code{setParam}
#' @param policy The optimal policy of the MDP
#' @param WatObs Given soil wate content data (if givenWatInfo==TRUE)
#' @param temData Temperature data
#' @param precData Precipitation data
#' @param iniTrueWat Initial soil-water content at day t=1
#'
#' @return A data frame containing the optimal actions and updated information obtained by Gaussian SSM
#' @export
optimalSearch<-function(param, policy, WatObs, temData, precData, iniTrueWat, givenWatInfo){

  soilWatObs<-c()
  if(!givenWatInfo){
    soilWatObs<-simWat(temData, precData, iniTrueWat, param)$soilWatObs
  }else{
    soilWatObs<-WatObs
  }

  paramGSSM<-setModel(param)

  meanPos<-DLMfilter(mod = paramGSSM, D = soilWatObs, Tem = temData, Pre = precData, W = paramGSSM$W, V = paramGSSM$V)$meanPos
  varPos<-DLMfilter(mod = paramGSSM, D = soilWatObs, Tem = temData, Pre = precData, W = paramGSSM$W, V = paramGSSM$V)$varPos
  sdPos<-sqrt(varPos)

  optAction<-c()
  weight<-c()
  operation<-c()
  dayLeft<-c()
  idxMW<-c()
  idxSW<-c()
  idxMP<-c()
  idxSP<-c()
  idxT<-c()
  idxP<-c()

  for(t in 1:param$tMax){
    if(t==1){
      opeNext=1
      dLNext=param$opD[opeNext]
      idxMW[t]<-findIndex(iniTrueWat,param$disAvgWat)
      idxSW[t]<-0
      idxMP[t]<-findIndex(param$gSSMm0,param$disMeanPos)
      idxSP[t]<-findIndex(sqrt(param$gSSMc0),param$disSdPos)
      idxT[t]<-findIndex(temData[t],param$disTem)
      idxP[t]<-findIndex(precData[t],param$disPre)
    }else{
      idxMW[t]<-findIndex(soilWatObs[t],param$disAvgWat)
      idxSW[t]<-0
      idxMP[t]<-findIndex(meanPos[t],param$disMeanPos)
      idxSP[t]<-findIndex(sdPos[t],param$disSdPos)
      idxT[t]<-findIndex(temData[t],param$disTem)
      idxP[t]<-findIndex(precData[t],param$disPre)
    }

    ope<-opeNext
    dL<-dLNext
    operation[t]<-ope
    dayLeft[t]<-dL
    optAction[t]<-subset(policy, day==t & op==ope & d==dL & iMW==idxMW[t] & iSW==idxSW[t] & iMP==idxMP[t] & iSP==idxSP[t] & iT==idxT[t] & iP==idxP[t]) ["optAction"]
    weight[t]<-as.numeric(subset(policy, day==t & op==ope & d==dL & iMW==idxMW[t] & iSW==idxSW[t] & iMP==idxMP[t] & iSP==idxSP[t] & iT==idxT[t] & iP==idxP[t]) ["weight"] )

    if(optAction[t]=="pos."){
      dLNext=dL
      opeNext=ope
    }
    if( ( (optAction[t]=="do.") || (optAction[t]=="doF.") ) && (dL>1)  ){
      dLNext=dL-1
      opeNext=ope
    }
    if( ( (optAction[t]=="do.") || (optAction[t]=="doF.") ) && (dL==1) && (ope<param$opNum) ){
      dLNext=param$opD[ope+1]
      opeNext=ope+1
    }
    if( ( (optAction[t]=="do.") || (optAction[t]=="doF.") ) && (dL==1) && (ope==param$opNum) ){
      tTerm<-t
      break
    }

  }

  dat<-data.table(t=1:(param$tMax))
  dat$MW<-soilWatObs[1:param$tMax]
#  dat$SW<-0
  dat$MP<-meanPos[1:param$tMax]
  dat$SP<-sdPos[1:param$tMax]
  dat$Tem<-temData[1:param$tMax]
  dat$Pre<-precData[1:param$tMax]

  dat$iMW<-idxMW[1:param$tMax]
#  dat$iSW<-0
  dat$iMP<-idxMP[1:param$tMax]
  dat$iSP<-idxSP[1:param$tMax]
  dat$iTem<-idxT[1:param$tMax]
  dat$iPre<-idxP[1:param$tMax]

  dat$operation<-operation[1:param$tMax]
  dat$dayLeft<-dayLeft[1:param$tMax]
  dat$optAction<-optAction[1:param$tMax]
  dat$weight<-weight[1:param$tMax]

  return(dat)
}

