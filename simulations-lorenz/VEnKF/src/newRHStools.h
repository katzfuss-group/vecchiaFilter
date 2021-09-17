#include "vectorTools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec getW( const arma::vec& X, const int&K );

arma::vec getXX(const arma::vec& X, const int& K, const arma::vec& W);

arma::vec RHS(const arma::vec& X, const int& K, const double& F);