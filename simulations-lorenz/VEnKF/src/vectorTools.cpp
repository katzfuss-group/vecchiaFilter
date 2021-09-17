#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



arma::uword mod(const int& i, const int& n){

  int mod = i % n;
  if(mod < 0) {
    mod += n;
  }
  return mod;
}


arma::vec mask(const arma::vec& X, arma::uword first, arma::uword last, bool truncate){
  
  arma::uword N = X.n_rows;
  arma::uvec maskV(N);
  if( first <= last ){
    maskV.fill(0);
    maskV.subvec( first, last ).fill(1);
  } else {
    maskV.fill(0);
    maskV.subvec( 0, last).fill(1);
    maskV.subvec( first, N-1 ).fill(1);
  }
  
  arma::vec masked = X % maskV;
  
  if(truncate){
    masked = masked( find( maskV ) );
  }
  
  return( masked );
}


