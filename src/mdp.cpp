#include "mdp.h"

// ===================================================

const double MDPV::ZERO = 1e-10; //4.656613e-10;   // equal 1/2147483647 (limit of int since trans pr is stored as int when build the hgf)

// ===================================================

MDPV::MDPV(const List paramModel){

  List rParam(paramModel);       // Get parameters in params
  opNum=as<int>(rParam["opNum"]);
  tMax=as<int>(rParam["tMax"]);
  opSeq=as<arma::vec>(rParam["opSeq"]);
  opE=as<arma::vec>(rParam["opE"]);
  opL=as<arma::vec>(rParam["opL"]);
  opD=as<arma::vec>(rParam["opD"]);
  opDelay=as<arma::vec>(rParam["opDelay"]);
  opFixCost=as<arma::vec>(rParam["opFixCost"]);
  watTh=as<arma::vec>(rParam["watTh"]);
  coefLoss=as<double>(rParam["coefLoss"]);
  priceYield=as<double>(rParam["priceYield"]);
  machCap=as<double>(rParam["machCap"]);
  yieldHa=as<double>(rParam["yieldHa"]);
  fieldArea=as<double>(rParam["fieldArea"]);
  coefTimeliness=as<double>(rParam["coefTimeliness"]);
  costSkip=as<double>(rParam["costSkip"]);

  watUpper=as<arma::vec>(rParam["watUpper"]);
  watLower=as<arma::vec>(rParam["watLower"]);
  stress=as<arma::vec>(rParam["stress"]);
  strength=as<arma::vec>(rParam["strength"]);
  weightCompletion=as<double>(rParam["weightCompletion"]);
  weightWorkable=as<double>(rParam["weightWorkable"]);
  weightTraffic=as<double>(rParam["weightTraffic"]);
  minOpt=as<int>(rParam["minOpt"]);
  maxOpt=as<int>(rParam["maxOpt"]);

  temMeanDry=as<double>(rParam["temMeanDry"]);
  temMeanWet=as<double>(rParam["temMeanWet"]);
  temVarDry=as<double>(rParam["temVarDry"]);
  temVarWet=as<double>(rParam["temVarWet"]);
  dryDayTh=as<double>(rParam["dryDayTh"]);
  precShape=as<double>(rParam["precShape"]);
  precScale=as<double>(rParam["precScale"]);
  prDryWet=as<double>(rParam["prDryWet"]);
  prWetWet=as<double>(rParam["prWetWet"]);

  hydroWatR=as<double>(rParam["hydroWatR"]);
  hydroWatS=as<double>(rParam["hydroWatS"]);
  hydroM=as<double>(rParam["hydroM"]);
  hydroKs=as<double>(rParam["hydroKs"]);
  hydroFi=as<double>(rParam["hydroFi"]);
  hydroLamba=as<double>(rParam["hydroLamba"]);
  hydroETa=as<double>(rParam["hydroETa"]);
  hydroETb=as<double>(rParam["hydroETb"]);
  hydroETx=as<double>(rParam["hydroETx"]);

  gSSMW=as<double>(rParam["gSSMW"]);
  gSSMV=as<double>(rParam["gSSMV"]);
  gSSMm0=as<double>(rParam["gSSMm0"]);
  gSSMc0=as<double>(rParam["gSSMc0"]);
  nGSSMm0=as<double>(rParam["nGSSMm0"]);
  nGSSMc0=as<double>(rParam["nGSSMc0"]);
  nGSSMK=as<double>(rParam["nGSSMK"]);

  check = as<bool>(rParam["check"]);
  rewRisk = as<bool>(rParam["rewRisk"]);

  dMP = as<arma::mat>(rParam["disMeanPos"]);
  dSP = as<arma::mat>(rParam["disSdPos"]);
  dMW = as<arma::mat>(rParam["disAvgWat"]);
  dSW = as<arma::mat>(rParam["disSdWat"]);
  dT = as<arma::mat>(rParam["disTem"]);
  dP = as<arma::mat>(rParam["disPre"]);

  sMP = as<arma::vec>(rParam["centerPointsMeanPos"]);
  sSP = as<arma::vec>(rParam["centerPointsSdPos"]);
  sMW = as<arma::vec>(rParam["centerPointsAvgWat"]);
  sSW = as<arma::vec>(rParam["centerPointsSdWat"]);
  sT = as<arma::vec>(rParam["centerPointsTem"]);
  sP = as<arma::vec>(rParam["centerPointsPre"]);

  sizeSMP = sMP.size();
  sizeSSP = sSP.size();
  sizeSMW = sMW.size();
  sizeSSW = sSW.size();
  sizeST = sT.size();
  sizeSP = sP.size();

  // matrices for filling the rewards and transition probabilities before running the HMDP:
  prMW = vector <vector<vector< vector< vector< vector<double> > > > > >(sizeSMW,
         vector<vector< vector< vector< vector<double> > > > >(sizeSMP,
         vector<vector< vector< vector<double> > > >(sizeSSP,
         vector<vector< vector<double> > >(sizeST,
         vector<vector<double> > (sizeSP,
         vector <double>(sizeSMW) ) ) ) ) ); //prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW]

  prMP = vector <vector<vector< vector< vector< vector<double> > > > > >(sizeSMW,
         vector<vector< vector< vector< vector<double> > > > >(sizeSMP,
         vector<vector< vector< vector<double> > > >(sizeSSP,
         vector<vector< vector<double> > >(sizeST,
         vector<vector<double> > (sizeSP,
         vector <double>(sizeSMP) ) ) ) ) ); //prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP]

  prSP = vector<vector< vector< vector< vector<double> > > > >(sizeSMW,
         vector<vector< vector< vector<double> > > >(sizeSSP,
         vector<vector< vector<double> > >(sizeST,
         vector<vector<double> > (sizeSP,
         vector <double>(sizeSSP) ) ) ) ); // prSP[iMWt][iSPt][iTt][iPt][iSP]


  prSW = vector < vector< vector<double> > > (tMax+1,
         vector< vector<double> >(sizeSSW,
         vector<double>(sizeSSW) ) ); //prSW[t][iSWt][iSW]

  prT =  vector <vector< vector<double> > > (sizeST,
         vector< vector<double> > (sizeSP,
         vector<double>(sizeST) ) ); //prT[iTt][iPt][iT]

  prP = vector< vector<double> > (sizeSP,
        vector<double>(sizeSP) ); //prP[iPt][iP]

  rewDo = vector <vector< vector<double> > >(opNum,
          vector< vector<double> >(sizeSMW,
          vector<double>(sizeSSW) ) ); //rewDo[op][iMWt][iSWt]


  int opDMax = arma::max(opD);

  mapL1Vector = vector< vector< vector< vector< vector< vector< vector< vector<int> > > > > > > >(opNum,
                vector< vector< vector< vector< vector< vector< vector<int> > > > > > >(opDMax+1,
                vector< vector< vector< vector< vector< vector<int> > > > > >(sizeSMW,
                vector< vector< vector< vector< vector<int> > > > >(sizeSSW,
                vector< vector< vector< vector<int> > > >(sizeSMP,
                vector< vector< vector<int> > >(sizeSSP,
                vector< vector<int> >(sizeST,
                vector <int>(sizeSP) ) ) ) ) ) ) ) ; //mapL1Vector[op][d][iMW][iSW][iMP][iSP][iT][iP];

  valueFun =    vector< vector< vector< vector< vector< vector< vector< vector< vector<double> > > > > > > > >(tMax+1,
                vector< vector< vector< vector< vector< vector< vector< vector<double> > > > > > > >(opNum,
                vector< vector< vector< vector< vector< vector< vector<double> > > > > > >(opDMax+1,
                vector< vector< vector< vector< vector< vector<double> > > > > >(sizeSMW,
                vector< vector< vector< vector< vector<double> > > > >(sizeSSW,
                vector< vector< vector< vector<double> > > >(sizeSMP,
                vector< vector< vector<double> > >(sizeSSP,
                vector< vector<double> >(sizeST,
                vector <double>(sizeSP) ) ) ) ) ) ) ) ) ; //valueFun[t][op][d][iMW][iSW][iMP][iSP][iT][iP];

  valFunDummy = vector<double>(tMax+1);

  optAction =   vector< vector< vector< vector< vector< vector< vector< vector< vector<string> > > > > > > > >(tMax+1,
                vector< vector< vector< vector< vector< vector< vector< vector<string> > > > > > > >(opNum,
                vector< vector< vector< vector< vector< vector< vector<string> > > > > > >(opDMax+1,
                vector< vector< vector< vector< vector< vector<string> > > > > >(sizeSMW,
                vector< vector< vector< vector< vector<string> > > > >(sizeSSW,
                vector< vector< vector< vector<string> > > >(sizeSMP,
                vector< vector< vector<string> > >(sizeSSP,
                vector< vector<string> >(sizeST,
                vector <string>(sizeSP) ) ) ) ) ) ) ) ) ; //optAction[t][op][d][iMW][iSW][iMP][iSP][iT][iP];

}

// ===================================================

void MDPV::Preprocess() {
  Rcout << "Build the HMDP ... \n\nStart preprocessing ...\n"<<endl;
  CalcTransPrMW();
  CalcTransPrSW();
  CalcTransPrMP();
  CalcTransPrSP();
  CalcTransPrT();
  CalcTransPrP();
  CalcRewaerdDo();
  Rcout << "... finished preprocessing.\n";
}



// ===================================================
SEXP MDPV::SolveMDP(){
  Preprocess();
  int t, op, iMW, iSW, iMP, iSP, iT, iP, d;
  double valueDo, valuePos;

  int counter=0;

  for(t=tMax; t>=1; --t){
    if(t==tMax){
      if (rewRisk) valFunDummy[tMax]=0; else valFunDummy[tMax]=priceYield*yieldHa*fieldArea;
      continue;
      }
    Rcout<<" day: "<<t<<endl;
    for(op=0; op<opNum; op++){
      if( (opE[op]>t) || (opL[op]<=t) ) continue;
      for(d=1; d<=opD[op]; d++){
        if(opD[op]-t+opE(op)>d) continue;
        if(opL[op]-t<d) continue;
        for(iMW=0; iMW<sizeSMW; iMW++){
          for(iSW=0; iSW<sizeSSW; iSW++){
            for(iMP=0; iMP<sizeSMP; iMP++){
              for(iSP=0; iSP<sizeSSP; iSP++){
                for(iT=0; iT<sizeST; iT++){
                  for(iP=0; iP<sizeSP; iP++){
                    if ( d<opL[op]-t ){
                      valuePos=WeightPos(op,d,iMW,iSW,iMP,iSP,iT,iP,t); counter = counter+1;
                      valueDo=WeightDo(op,d,iMW,iSW,iMP,iSP,iT,iP,t); counter = counter+1;
                      if(valueDo>valuePos){
                        valueFun[t][op][d][iMW][iSW][iMP][iSP][iT][iP]=valueDo; optAction[t][op][d][iMW][iSW][iMP][iSP][iT][iP]="do.";
                      }else{
                        valueFun[t][op][d][iMW][iSW][iMP][iSP][iT][iP]=valuePos; optAction[t][op][d][iMW][iSW][iMP][iSP][iT][iP]="pos.";
                      }
                    }
                    if( d==(opL[op]-t) ){
                      valueDo=WeightDo(op,d,iMW,iSW,iMP,iSP,iT,iP,t); counter = counter+1;
                      valueFun[t][op][d][iMW][iSW][iMP][iSP][iT][iP]=valueDo; optAction[t][op][d][iMW][iSW][iMP][iSP][iT][iP]="doF.";
                    }
                    valFunDummy[t]=0+valFunDummy[t+1]; //IS IT TRUE?

                  }
                }
              }
            }
          }
        }
      }
    }
  }
  Rcout<<" Number of actions: "<< counter << endl;
  totalRew=weightIni();
  printPolicy();
  return( wrap( List::create(Named("totalRew") = totalRew) ) );
  // return( wrap( List::create(Named("weights") = valueFun, Named("optAction") = optAction, Named("totalRew") = totalRew) ) );

}


// ===================================================

double MDPV::WeightPos(int & opt, int & dt, int & iMWt, int & iSWt, int & iMPt, int & iSPt, int & iTt, int & iPt, int & t) {
  double pr4, prS, reward;
  double weightFu=0;
  int op,d,iMW,iSW,iMP,iSP,iT,iP, tN;
  op=opt;
  d=dt;
  tN=t+1;

  pr.clear();
  for(iMW=0; iMW<sizeSMW; iMW++){
    for(iSW=0; iSW<sizeSSW; iSW++){
      for(iMP=0; iMP<sizeSMP; iMP++){
        for(iSP=0; iSP<sizeSSP; iSP++){
          for(iT=0; iT<sizeST; iT++){
            for(iP=0; iP<sizeSP; iP++){
              prS = prSP[iMWt][iSPt][iTt][iPt][iSP];
              pr4 = prS*exp(prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW] + prSW[t][iSWt][iSW] + prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP]
                              + prT[iTt][iPt][iT] + prP[iPt][iP]);
              if (pr4>ZERO) {
                weightFu = weightFu + pr4*valueFun[tN][op][d][iMW][iSW][iMP][iSP][iT][iP];
                pr.push_back(pr4);
              }
            }
          }
        }
      }
    }
  }

  if(rewRisk) reward=0; else  reward=-coefTimeliness*priceYield*yieldHa*fieldArea;

  if (check) {
    arma::vec tmp(pr);
    if (!Equal(sum(tmp),1,1e-8)) {
      Rcout << "Warning sum pr!=1 in WeightsTransPrPos - diff = " << 1-sum(tmp) << " op = " << op << " action = pos. " << " index:" << endl; //vec2String<int>(index) << " pr:" << vec2String<flt>(pr) << endl;
    }
  }
  return(reward+weightFu);
}

// ===================================================

double MDPV::WeightDo(int & opt, int & dt, int & iMWt, int & iSWt, int & iMPt, int & iSPt, int & iTt, int & iPt, int & t) {
  double pr4, prS, reward;
  double weightFu=0;
  double completionCri=0;
  int op,d,iMW,iSW,iMP,iSP,iT,iP, tN;

  pr.clear();
  if( (dt>1) ){
    d=dt-1;
    op=opt;
    tN=t+1;
    for(iMW=0; iMW<sizeSMW; iMW++){
      for(iSW=0; iSW<sizeSSW; iSW++){
        for(iMP=0; iMP<sizeSMP; iMP++){
          for(iSP=0; iSP<sizeSSP; iSP++){
            for(iT=0; iT<sizeST; iT++){
              for(iP=0; iP<sizeSP; iP++){
                prS = prSP[iMWt][iSPt][iTt][iPt][iSP];
                pr4 = prS*exp(prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW] + prSW[t][iSWt][iSW] + prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP]
                                + prT[iTt][iPt][iT] + prP[iPt][iP]);
                if (pr4>ZERO) {
                  weightFu = weightFu + pr4*valueFun[tN][op][d][iMW][iSW][iMP][iSP][iT][iP];
                  pr.push_back(pr4);
                }
              }
            }
          }
        }
      }
    }

    if(rewRisk) reward = rewDo[opt][iMWt][iSWt]; else reward=rewDo[opt][iMWt][iSWt];
  }

  if( (dt==1) & (opt<(opNum-1)) ){
    d=opD[opt+1];
    op=opt+1;
    tN=t+1;
    for(iMW=0; iMW<sizeSMW; iMW++){
      for(iSW=0; iSW<sizeSSW; iSW++){
        for(iMP=0; iMP<sizeSMP; iMP++){
          for(iSP=0; iSP<sizeSSP; iSP++){
            for(iT=0; iT<sizeST; iT++){
              for(iP=0; iP<sizeSP; iP++){
                prS = prSP[iMWt][iSPt][iTt][iPt][iSP];
                pr4 = prS*exp(prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW] + prSW[t][iSWt][iSW] + prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP]
                                + prT[iTt][iPt][iT] + prP[iPt][iP]);
                if (pr4>ZERO) {
                  weightFu = weightFu + pr4*valueFun[tN][op][d][iMW][iSW][iMP][iSP][iT][iP];
                  pr.push_back(pr4);
                }
              }
            }
          }
        }
      }
    }
    if(rewRisk) reward = rewDo[opt][iMWt][iSWt]; else reward=rewDo[opt][iMWt][iSWt];
  }

  if( (dt==1) & (opt==(opNum-1)) ){
    pr4=1;
    pr.push_back(pr4);
    tN=t+1;
    weightFu = weightFu + pr4*valFunDummy[tN];
    if( (tN>=minOpt) & (tN<=maxOpt) ) completionCri = 1;
    if( tN<minOpt ) completionCri = (double)(minOpt-tN)/(double)(minOpt);
    if( tN>maxOpt ) completionCri = (double)(tN-maxOpt)/(double)(tN);

    if(rewRisk) reward = ( rewDo[opt][iMWt][iSWt] +  weightCompletion*(completionCri) ); else reward=rewDo[opt][iMWt][iSWt];
  }

  if (check) {
    arma::vec tmp(pr);
    if (!Equal(sum(tmp),1,1e-8)) {
      Rcout << "Warning sum pr!=1 in WeightsTransPrDo - diff = " << 1-sum(tmp) << " op = " << op << " action = Do. " << endl; // " index:" << vec2String<int>(index) << " pr:" << vec2String<flt>(pr) << endl;
    }
  }

  return(reward+weightFu);
}

// ===================================================

double MDPV::weightIni() {
  // double pr4, prS, reward;
  // double weightFu=0;
  // int op,d,iMW,iSW,iMP,iSP,iT,iP, id, tN;
  // int iMWt,iSWt,iMPt,iSPt,iTt,iPt;
  // iPt= findIndex(precShape*precScale,dP);
  // if(precShape*precScale<0.25) iTt= findIndex(temMeanDry,dT); else iTt= findIndex(temMeanWet,dT);
  // op=0;
  // d=opD[0];
  // tN=1;
  //
  // pr.clear();
  // for(iMW=0; iMW<sizeSMW; iMW++){
  //   for(iSW=0; iSW<sizeSSW; iSW++){
  //     for(iMP=0; iMP<sizeSMP; iMP++){
  //       for(iSP=0; iSP<sizeSSP; iSP++){
  //         for(iT=0; iT<sizeST; iT++){
  //           for(iP=0; iP<sizeSP; iP++){
  //             prS = prSP[iMWt][iSPt][iTt][iPt][iSP];
  //             pr4 = prS*exp(prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW] + prSW[t][iSWt][iSW] + prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP]
  //                             + prT[iTt][iPt][iT] + prP[iPt][iP]);
  //             if (pr4>ZERO) {
  //               weightFu = weightFu + pr4*valueFun[tN][op][d][iMW][iSW][iMP][iSP][iT][iP];
  //               pr.push_back(pr4);
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  // reward=0
  //
  // if (check) {
  //   arma::vec tmp(pr);
  //   if (!Equal(sum(tmp),1,1e-8)) {
  //     Rcout << "Warning sum pr!=1 in WeightsTransPrPos - diff = " << 1-sum(tmp) << " op = " << op << " action = pos. " << " index:" << endl; //vec2String<int>(index) << " pr:" << vec2String<flt>(pr) << endl;
  //   }
  // }
  // return(reward+weightFu);
  return(0);
}


// ===================================================

void MDPV::CalcRewaerdDo(){
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int iMW,iSW,op;
  double workCri,trafiCriteria;

  for(op=0; op<opNum; op++){
    for(iMW=0;iMW<sizeSMW;iMW++){
      for(iSW=0;iSW<sizeSSW;iSW++){
        if(rewRisk){
          if( (dMW(iMW,0)>=watLower[op] ) & (dMW(iMW,0)<=watUpper[op] ) ) workCri=1;
          if( dMW(iMW,0)<watLower[op] ) workCri=(double)(watLower[op]-dMW(iMW,0))/(double)(watLower[op]);
          if( dMW(iMW,0)>watUpper[op] ) workCri=(double)(dMW(iMW,0)-watUpper[op])/(double)(dMW(iMW,0));

          if(strength[iMW]>=stress[iMW]) trafiCriteria=1;
          if(strength[iMW]<stress[iMW]) trafiCriteria=(double)(stress[iMW]-strength[iMW])/(double)(stress[iMW]);

          rewDo[op][iMW][iSW] = weightWorkable*workCri + weightTraffic*trafiCriteria;
        }else{
          rewDo[op][iMW][iSW]= -coefLoss*priceYield*yieldHa*machCap*(1- R::pnorm(watTh[op],dMW(iMW,0),dSW(iSW,0),1,0) );
        }
      }
    }
  }
}

// ===================================================
void MDPV::CalcTransPrMW(){   //prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW]
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int iMWt, iMPt,iSPt,iTt,iPt,iMW;
  double lower, upper, mt, ct, ft,rt,qt;

  for(iMWt=0; iMWt<sizeSMW; iMWt++){
    for(iMPt=0;iMPt<sizeSMP;iMPt++){
      for(iSPt=0;iSPt<sizeSSP;iSPt++){
        for(iTt=0;iTt<sizeST;iTt++){
          for(iPt=0;iPt<sizeSP;iPt++){
            ft = Hydro(dMW(iMWt,0),dT(iTt,0),dP(iPt,0));
            rt= pow(dSP(iSPt,0),2) + gSSMW;
            qt=pow(ft,2)*rt + gSSMV;
            mt=ft*dMP(iMPt,0); ct=qt;
            for(iMW=0; iMW<sizeSMW; iMW++){
              lower= dMW(iMW,1); upper=dMW(iMW,2);
              prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW] = log( R::pnorm(upper,mt,sqrt(ct),1,0) - R::pnorm(lower,mt,sqrt(ct),1,0)  );
            }
          }
        }
      }
    }
  }
}

// ===================================================

void MDPV::CalcTransPrSW(){  //prSW[t][iSWt][iSW]
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int t,iSWt,iSW;

  for(t=1; t<tMax; t++){
    for(iSWt=0; iSWt<sizeSSW; iSWt++){
      for(iSW=0; iSW<sizeSSW; iSW++){
        prSW[t][iSWt][iSW] = log(1);//log(PrNGSSM( t, pow(dSW(iSW,1),2), pow(dSW(iSW,2),2), pow(dSW(iSWt,0),2) ) );
      }
    }
  }
}


// ===================================================

void MDPV::CalcTransPrMP(){   //prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP]
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int iMWt, iMPt,iSPt,iTt,iPt,iMP;
  double lower, upper, mt, ct, ft, rt, qt;

  for(iMWt=0; iMWt<sizeSMW; iMWt++){
    for(iMPt=0;iMPt<sizeSMP;iMPt++){
      for(iSPt=0;iSPt<sizeSSP;iSPt++){
        for(iTt=0;iTt<sizeST;iTt++){
          for(iPt=0;iPt<sizeSP;iPt++){
            ft = Hydro(dMW(iMWt,0),dT(iTt,0),dP(iPt,0));
            rt= pow(dSP(iSPt,0),2) + gSSMW;
            qt=pow(ft,2)*rt + gSSMV;
            mt=dMP(iMPt,0); ct= pow(rt,2)*pow(ft,2)/qt;
            for(iMP=0; iMP<sizeSMP; iMP++){
              lower= dMP(iMP,1); upper=dMP(iMP,2);
              prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP] = log( R::pnorm(upper,mt,sqrt(ct),1,0) - R::pnorm(lower,mt,sqrt(ct),1,0)  );
            }
          }
        }
      }
    }
  }
}
// ===================================================

void MDPV::CalcTransPrSP(){ //prSP[iMWt][iSPt][iTt][iPt][iSP]
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int iMWt,iSPt,iTt,iPt,iSP;
  double ft, ct, cN;

  for(iMWt=0;iMWt<sizeSMW;iMWt++){
    for(iSPt=0;iSPt<sizeSSP;iSPt++){
      for(iTt=0;iTt<sizeST;iTt++){
        for(iPt=0;iPt<sizeSP;iPt++){
          ft=Hydro(dMW(iMWt,0),dT(iTt,0),dP(iPt,0));
          ct=pow(dSP(iSPt,0),2);
          cN= ( (ct+gSSMW)*gSSMV )/( pow(ft,2)*(ct+gSSMW) + gSSMV);
          for(iSP=0;iSP<sizeSSP;iSP++){
            if( ( cN>pow(dSP(iSP,1),2) ) & ( cN<=pow(dSP(iSP,2),2) ) ){
              prSP[iMWt][iSPt][iTt][iPt][iSP]=1;
            }else{
              prSP[iMWt][iSPt][iTt][iPt][iSP]=0;
            }
          }
        }
      }
    }
  }
}

// ===================================================

void MDPV::CalcTransPrT(){ //prT[iTt][iPt][iT]
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int iTt,iT,iPt;
  double lower, upper, mt, ct;

  for(iTt=0; iTt<sizeST; iTt++){
    for(iPt=0; iPt<sizeSP; iPt++){
      if(dP(iPt,0)<=dryDayTh){
        mt=temMeanDry; ct=temVarDry;
      }else{
        mt=temMeanWet; ct=temVarWet;
      }
      for(iT=0; iT<sizeST; iT++){
        lower=dT(iT,1); upper=dT(iT,2);
        prT[iTt][iPt][iT] = log( R::pnorm(upper,mt,sqrt(ct),1,0) - R::pnorm(lower,mt,sqrt(ct),1,0) );
      }
    }
  }
}
// ===================================================

void MDPV::CalcTransPrP(){ //prP[iPt][iP]
  cpuTime.Reset(0); cpuTime.StartTime(0);
  int iPt,iP;
  double lower, upper; //, lowert, uppert;

  for(iPt=0; iPt<sizeSP; iPt++){
    //uppert=dP(iPt,1); lowert=dP(iPt,2);
    for(iP=0; iP<sizeSP; iP++){
      lower=dP(iP,1); upper=dP(iP,2);
      //if( (uppert<=dryDayTh) & (upper<=dryDayTh) ) prP[iPt][iP]= log(1-prDryWet);
      //if( (lowert>dryDayTh) & (upper<=dryDayTh) ) prP[iPt][iP]= log(1-prWetWet);
      //if( (uppert<=dryDayTh) & (lower>dryDayTh) ) prP[iPt][iP]= log(prDryWet*( R::pgamma(upper,precShape,precScale,1,0) - R::pgamma(lower,precShape,precScale,1,0)) );
      //if( (lowert>dryDayTh) & (lower>dryDayTh) ) prP[iPt][iP]= log(prWetWet*( R::pgamma(upper,precShape,precScale,1,0) - R::pgamma(lower,precShape,precScale,1,0) ) );

      // if( (dP(iPt,0)<=dryDayTh) & (dP(iP,0)<=dryDayTh) ) prP[iPt][iP]= log(1-prDryWet);
      // if( (dP(iPt,0)>dryDayTh) & (dP(iP,0)<=dryDayTh) ) prP[iPt][iP]= log(1-prWetWet);
      // if( (dP(iPt,0)<=dryDayTh) & (dP(iP,0)>dryDayTh) ) prP[iPt][iP]= log(prDryWet*( R::pgamma(upper,precShape,precScale,1,0) - R::pgamma(lower,precShape,precScale,1,0)) );
      // if( (dP(iPt,0)>dryDayTh) & (dP(iP,0)>dryDayTh) ) prP[iPt][iP]= log(prWetWet*( R::pgamma(upper,precShape,precScale,1,0) - R::pgamma(lower,precShape,precScale,1,0) ) );

      if( (dP(iPt,0)<=dryDayTh) & (dP(iP,0)<=dryDayTh) ) prP[iPt][iP]= log(1-prDryWet);
      if( (dP(iPt,0)>dryDayTh) & (dP(iP,0)<=dryDayTh) ) prP[iPt][iP]= log(1-prWetWet);
      if( (dP(iPt,0)<=dryDayTh) & (dP(iP,0)>dryDayTh)  ) { if(iP==1)lower=0; prP[iPt][iP]= log(prDryWet*( R::pgamma(upper,precShape,precScale,1,0) - R::pgamma(lower,precShape,precScale,1,0)) ); }
      if( (dP(iPt,0)>dryDayTh) & (dP(iP,0)>dryDayTh) ) { if(iP==1)lower=0; prP[iPt][iP]= log(prWetWet*( R::pgamma(upper,precShape,precScale,1,0) - R::pgamma(lower,precShape,precScale,1,0) ) ); }
    }
  }
}

// ===================================================

double MDPV::Hydro(double & Wt, double & Tt, double Pt){ // ssm.hydro(dMW(iMWt,0),dT(iTt,0),dP(iPt,0));

  double f, e, g, ET;

  ET= hydroETa + hydroETb*hydroETx*(0.46*Tt + 8.13);
  //f= hydroKs*( 1-(hydroFi*(hydroWatS-Wt)/hydroF) );
  f=Pt*(1- pow((Wt-hydroWatR)/(hydroWatS-hydroWatR),hydroM) );
  g= hydroKs*pow((Wt-hydroWatR)/(hydroWatS-hydroWatR),3+2/hydroLamba);
  e= ET*(Wt-hydroWatR)/(hydroWatS-hydroWatR);

  return(Wt+f-e-g);
}


// ===================================================
double MDPV::PrNGSSM(int t, double lower, double upper, double var) { //SOLVED[Reza] : Based on the formulatiom for this probability in the paper, I changed "n" to "nf" (nf is the sample size).

  double probSd, xUpper, xLower;
  double a, s, alpha, gamma, beta;
  double G=1;

  double oShape = (double) (nGSSMK-1)/(2);   //shape parameter of observation distribution
  double iShape = (double) (nGSSMK-1)/(2);//(double) (numSample-3)/(numSample-5); //("shape parameter of prior at t=1"): c_1 in the paper

  a = (double) (G * var * ( iShape + oShape*t ) ) / ( iShape + oShape*(t+1) ) ; // location
  s = a; // scale
  alpha = oShape; // shape 1
  gamma = iShape + oShape*t + 1; // shape 2
  beta =1 ;  // Weibul parameter
  xUpper = (double) (1)/( 1 + pow ( (double)(upper-a)/(s),-beta ) );
  xLower = (double) (1)/( 1 + pow ( (double)(lower-a)/(s),-beta ) );

  probSd= R::pbeta(xUpper, alpha, gamma,1, 0) - R::pbeta(xLower, alpha, gamma,1, 0);
  if( ( R::pbeta(xUpper, alpha, gamma,1, 0) - R::pbeta(xLower, alpha, gamma,1, 0) )<0 ) DBG4("error_minus"<<endl)
    //if(probSd!=probSd)  DBG4(endl << " Error" << " lower=" << lower << " upper=" << upper << " centerp=" << var <<" t: "<<t<< endl)
    return (probSd);  // Rf_pgamma(q, shape, scale, lower.tail, log.p)
}


// ===================================================

int MDPV::countStatesMDP(){
  int x, t, iMW, iSW, iMP, iSP, iT, iP, op, d;
  x=0;

  for(t=tMax; t>=1; --t){
    if(t == tMax)  continue;
    for(op=0; op<opNum; op++){
      if( (opE[op]>t) || (opL[op]<=t) ) continue;
      for(d=1; d<=opD[op]; d++){
        if(opD[op]-t+opE(op)>d) continue;
        if(opL[op]-t<d) continue;
        for(iMW=0; iMW<sizeSMW; iMW++){
          for(iSW=0; iSW<sizeSSW; iSW++){
            for(iMP=0; iMP<sizeSMP; iMP++){
              for(iSP=0; iSP<sizeSSP; iSP++){
                for(iT=0; iT<sizeST; iT++){
                  for(iP=0; iP<sizeSP; iP++){
                    x++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(x+tMax);
}

// ===================================================


void MDPV::printPolicy(){
  int t, op, iMW, iSW, iMP, iSP, iT, iP, d;

  //Store the resalts in the csv files:
  ofstream  myFile;
  //myFile.open("C:\\Academic_Postdoc\\Codes\\hmdpTillage\\policyMDP.csv", ios::trunc);
  myFile.open("policyMDP.csv", ios::trunc);
  myFile << "statLbl" << ";" << "day" << ";" << "op" << ";" << "d" << ";" << "iMW" << ";" << "iSW" << ";" << "iMP" << ";" << "iSP" << ";" << "iT" << ";" << "iP" << ";" << "optAction" << ";" << "weight" <<endl;


  for(t=tMax; t>=1; --t){
    if(t==tMax){valFunDummy[tMax]=priceYield*yieldHa*fieldArea; continue;}
    for(op=0; op<opNum; op++){
      if( (opE[op]>t) || (opL[op]<=t) ) continue;
      for(d=1; d<=opD[op]; d++){
        if(opD[op]-t+opE(op)>d) continue;
        if(opL[op]-t<d) continue;
        for(iMW=0; iMW<sizeSMW; iMW++){
          for(iSW=0; iSW<sizeSSW; iSW++){
            for(iMP=0; iMP<sizeSMP; iMP++){
              for(iSP=0; iSP<sizeSSP; iSP++){
                for(iT=0; iT<sizeST; iT++){
                  for(iP=0; iP<sizeSP; iP++){
                    label = getLabel(op,d,iMW,iSW,iMP,iSP,iT,iP,t);
                    myFile << label << ";" << t << ";" << op + 1 << ";" << d << ";" << iMW << ";" <<
                              iSW << ";" << iMP << ";" << iSP << ";" << iT << ";" << iP << ";" <<
                              optAction[t][op][d][iMW][iSW][iMP][iSP][iT][iP] << ";" <<
                              valueFun[t][op][d][iMW][iSW][iMP][iSP][iT][iP] <<endl;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  myFile.close();
}

// ===================================================

int MDPV::findIndex(double st, arma::mat dis){
  for(int i=0; i<dis.n_rows; i++){
    if( ( st>=dis(i,1) ) & ( st<dis(i,2) )  )
      return(i);
  }
  cout<<"error in index: "<< st << " dis " << dis << endl;
  return(-1);
}


