#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::uword mod(const int&, const int&);

arma::vec mask(const arma::vec&, arma::uword, arma::uword, bool);