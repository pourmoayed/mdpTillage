
#ifndef MDPV_HPP
#define MDPV_HPP

#include "RcppArmadillo.h"    // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "binaryMDPWriter.h"
#include "time.h"

using namespace Rcpp;
using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do

// [[Rcpp::depends(RcppArmadillo)]]

// ===================================================

/**
* Class for soving an MDP model using value iteration algorithm for scheduling tillage operations.
*
* @author Reza Pourmoayed
*/
class MDPV
{
  public:  // methods

    /** Constructor. Store the parameters.
    *
    * @param paramModel A list for model parameters related to HMDP and SSMs created using \code{setParam} in R.
    */
    MDPV(const List paramModel);


    /** Build the HMDP (to binary files). Use "shared linking".
    *
    *  Build a 3 level HMDP saved in binary files by using the
    *  binaryMDPWriter (c++ version).
    *
    *  @return Build log (string)
    *  @export
    */
    SEXP SolveMDP();


    /** Count the number of states in the HMDP */
    int countStatesMDP();


private:

  /** Calculate and fill arrays with rewards and trans pr. */
  void Preprocess();



  /** Calculate the value function for action "pos." related to postpone tillage operation.
  *
  * @param op Tillage operation under consideration
  * @param dt Index of state for remaining days for finishing operation op
  * @param iMWt Index of state for estimated mean of soil water content
  * @param iSWt Index of state for estimated standard deviation of soil water content
  * @param iMPt Index of state for estimated posterior mean of latent variable in Gaussian SSM (error factor)
  * @param iSPt Index of state for estimated posterior sandard deviation of latent variable in Gaussian SSM (error factor)
  * @param iTt Index of state for weather forecast regarding air temprature.
  * @param iPt Index of state for weather forecast regarding precipitation.
  * @param t Current day.
  *
  */
  double WeightPos(int & op, int & dt, int & iMWt, int & iSWt, int & iMPt, int & iSPt, int & iTt, int & iPt, int & t);


  /** Calculate the value function for action "do." related to performing a tillage operation.
  *
  * @param op Tillage operation under consideration
  * @param dt Index of state for remaining days for finishing operation op
  * @param iMWt Index of state for estimated mean of soil water content
  * @param iSWt Index of state for estimated standard deviation of soil water content
  * @param iMPt Index of state for estimated posterior mean of latent variable in Gaussian SSM (error factor)
  * @param iSPt Index of state for estimated posterior sandard deviation of latent variable in Gaussian SSM (error factor)
  * @param iTt Index of state for weather forecast regarding air temprature.
  * @param iPt Index of state for weather forecast regarding precipitation.
  * @param t Current day.
  *
  */
  double WeightDo(int & op, int & dt, int & iMWt, int & iSWt, int & iMPt, int & iSPt, int & iTt, int & iPt, int & t);


  /** Calculate the total reward of the model (value function at the stage 0)
   */
  double weightIni();


  /** Calculate the reward values under action Do.
  *
  *  Values are stored in the vector \var(rewDo[op][iMW][iSW]).
  */
  void CalcRewaerdDo();

  /** Calculate the transition probability values for estimated mean of soil water content.
  *
  *  Values are stored in the vector \var(prMW[iMWt][iMPt][iSPt][iTt][iPt][iMW]).
  */
  void CalcTransPrMW();


  /** Calculate the transition probability values for estimated standard deviation of soil water content.
  *
  *  Values are stored in the vector \var(prSW[t][iSWt][iSW]).
  */
  void CalcTransPrSW();


  /** Calculate the transition probability values for posterior mean of latent variable (error factor) in Gaussian SSM.
  *
  *  Values are stored in the vector \var(prMP[iMWt][iMPt][iSPt][iTt][iPt][iMP]).
  */
  void CalcTransPrMP();


  /** Calculate the transition probability values for posterior standard deviation of latent variable (error factor) in Gaussian SSM.
  *
  *  Values are stored in the vector \var(prSP[iMWt][iSPt][iTt][iPt][iSP]).
  */
  void CalcTransPrSP();


  /** Calculate the transition probability values for weather information regarding air temprature.
  *
  *  Values are stored in the vector \var(prT[iTt][iPt][iT]).
  */
  void CalcTransPrT();


  /** Calculate the transition probability values for weather information regarding precipitation.
  *
  *  Values are stored in the vector \var(prP[iPt][iP]).
  */
  void CalcTransPrP();

  /** Print the optimal policy with optimal valur functions in a csv file.
   *
   */
  void printPolicy();


  /** Calculate the future soil water content based on a rainfall-runoff model given in \url(http://onlinelibrary.wiley.com/doi/10.1002/hyp.6629/abstract).
  *
  * @param Wt The current Soil water content (Volumetric measure)
  * @param Tt Prediction of average air temprature for the next day.
  * @param Pt Prediction of total precipitation for the next day.
  *
  * @return A predition of soil water content at the next day.
  */
  double Hydro(double & Wt, double & Tt, double Pt);


  /** Calculate the transition probability for variance component (posterior mean of latent variable) in a non-Gaussian SSM baesd on Theorem 4 in \url(http://www.sciencedirect.com/science/article/pii/S0377221715008802).
  *
  * @param t current day.
  * @param lower Lower limit of a variance component in the next day.
  * @param upper Upper limit of a variance component in the next day.
  * @param var Center point of the variance component in the current day.
  *
  * @return Transition probability for the variance component at the next day given variance var.
  */
  double PrNGSSM(int t, double lower, double upper, double var);

  /** Convert integers into a string. */
  string getLabel(const int & a, const int & b, const int & c, const int & d, const int & e, const int & f, const int & g, const int & h) {
    std::ostringstream s;
    s << "(" << a << "," << b << "," << c << "," << d << "," << e << "," << f << "," << g << "," << h << ")";
    return s.str();
  }

  /** Convert integers into a string. */
  string getLabel(const int & a, const int & b, const int & c, const int & d, const int & e, const int & f, const int & g, const int & h, const int & k) {
    std::ostringstream s;
    s << "(" << a << "," << b << "," << c << "," << d << "," << e << "," << f << "," << g << "," << h << "," << k << ")";
    return s.str();
  }

  /** Find the index of a given state based on the discritization matrix .
   *
   * @param st A value in a given discretized matrix
   * @param dis A matrix containing discretized values of a continuous state variable.
   *
   * @return The index number of a specefic interval in dis that includes st.
   */
  int findIndex(double st, arma::mat dis);

  //----------------------------------------------------------------------------------------------------------------------------------

  private:   // variables

    static const double ZERO;  // trans pr below are considered as ZERO

    arma::vec opSeq;
    arma::vec opE;
    arma::vec opL;
    arma::vec opD;
    arma::vec opDelay;
    arma::vec opFixCost;
    arma::vec watTh;

    arma::vec watUpper;
    arma::vec watLower;
    arma::vec stress;
    arma::vec strength;
    double weightCompletion;
    double weightWorkable;
    double weightTraffic;
    int minOpt;
    int maxOpt;

    int opNum;
    int tMax;
    double coefLoss;
    double priceYield;
    double machCap;
    double yieldHa;
    double fieldArea;
    double coefTimeliness;
    double costSkip;

    double temMeanDry;
    double temMeanWet;
    double temVarDry;
    double temVarWet;
    double dryDayTh;
    double precShape;
    double precScale;
    double prDryWet;
    double prWetWet;

    double hydroWatR;
    double hydroWatS;
    double hydroM;
    double hydroKs;
    double hydroFi;
    double hydroLamba;
    double hydroETa;
    double hydroETb;
    double hydroETx;

    double gSSMW;
    double gSSMV;
    double gSSMm0;
    double gSSMc0;
    double nGSSMm0;
    double nGSSMc0;
    double nGSSMK;

    bool check;
    bool rewRisk;

    arma::mat dMP;
    arma::mat dSP;
    arma::mat dMW;
    arma::mat dSW;
    arma::mat dT;
    arma::mat dP;

    arma::vec sMP;
    arma::vec sSP;
    arma::vec sMW;
    arma::vec sSW;
    arma::vec sT;
    arma::vec sP;

    int sizeSMP;
    int sizeSSP;
    int sizeSMW;
    int sizeSSW;
    int sizeST;
    int sizeSP;

    double totalRew;

    vector<flt> pr;


    vector <vector<vector< vector< vector< vector<double> > > > > > prMW;
    vector <vector<vector< vector< vector< vector<double> > > > > > prMP;
    vector<vector< vector< vector< vector<double> > > > > prSP;
    vector< vector< vector<double> > > prSW;
    vector <vector< vector<double> > > prT;
    vector< vector<double> > prP;
    vector <vector< vector<double> > > rewDo;
    vector< vector< vector< vector< vector< vector< vector< vector<int> > > > > > > > mapL1Vector;
    vector< vector< vector< vector< vector< vector< vector< vector< vector<double> > > > > > > > > valueFun;
    vector< vector< vector< vector< vector< vector< vector< vector< vector<string> > > > > > > > > optAction;
    vector<double> valFunDummy;

    TimeMan cpuTime;
    string label;
};


#endif
