#include "vectorTools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double Wn_Even_Cpp(const arma::vec& XX, const int& n, const int& k, const int& N);

double Wn_Odd_Cpp(const arma::vec& XX, const int& n, const int& k, const int& N);

double XX_Kn_Even_Cpp(const arma::vec& XX, const int& n, const int& k);

double XX_Kn_Odd_Cpp(const arma::vec& XX, const int& n, const int& k);

double oldRHS(const arma::vec& X, const int& j, const int& K, const double& F);