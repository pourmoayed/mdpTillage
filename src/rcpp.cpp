#include "mdp.h"

using namespace Rcpp;
using namespace std;

//' Solving the MDP using value iteration algorithm.
//'
//' @param paramModel parameters a list created using \code{\link{setParameters}}.
//'
//' @return A list
//' @export
// [[Rcpp::export]]
SEXP SolveMDPModel(const List paramModel) {
   MDPV Model(paramModel);
   Rcout << "Total number of states: " << Model.countStatesMDP() << endl;
   return( Model.SolveMDP() );
   //return(wrap(0));
}













