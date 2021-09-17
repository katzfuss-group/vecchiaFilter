#include "vectorTools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec getW( const arma::vec& X, const int&K ){
  
  arma::uword N = X.n_rows;
  arma::vec W(N);

  arma::uword J = floor(K/2);
  arma::uword n = 0;
  
  while( n < N ){
    
    arma::uword first = mod(n - J, N);
    arma::uword last = mod(n + J, N);
    arma::vec masked = mask(X, first, last, false);
    masked /= K;
    
    if( K % 2 == 0 ){
      masked(first) *= 0.5;
      masked(last) *= 0.5;  
    }
    W(n) = sum(masked);
    
    //Rcout << "first=" << first << std::endl;
    //Rcout << "last=" << last << std::endl;
    //Rcout << "mask: " << mask(X, first, last).t() << std::endl;
    //Rcout << "masked: " << masked.t() << std::endl;
    
    n++;  
    
  }
  
  return(W);
  
}


arma::vec getXX(const arma::vec& X, const int& K, const arma::vec& W){
  
  arma::uword N = X.n_rows;

  arma::vec XX = -shift(W, -2*K) % shift(W, -K);
  
  arma::uword J = floor(K/2);

  arma::vec dotProds = arma::zeros<arma::vec>(N);
  
  for( arma::uword n = 0; n < N; n++ ){
    
    arma::uword firstW = mod(n - K - J, N);
    arma::uword lastW = mod(n - K + J, N);
    arma::vec Wvec = mask(W, firstW, lastW, true);    
    arma::uword firstX = mod(n + K - J, N);
    arma::uword lastX = mod(n + K + J, N);
    arma::vec Xvec = mask(X, firstX, lastX, true);
    
    if(K % 2==0){
      Wvec.head(1) /= sqrt(2);
      Wvec.tail(1) /= sqrt(2);  
      
      Xvec.head(1) /= sqrt(2);
      Xvec.tail(1) /= sqrt(2);
    }

    //Rcout << "masked W:" << Wvec.t() << std::endl;
    //Rcout << "masked X:" << Xvec.t() << std::endl;

    arma::vec dp = Xvec.t() * Wvec/K;
    //Rcout << "dot product: " << dp << std::endl;
    dotProds(n) = dp(0,0);
  }
  
  return( XX + dotProds );
  
}


//' right hand side of the Lorenz model in vector form
//'
//' @param X current state
//' @param K parameter from the Lorenz 04 model
//' @param F parameter from the Lorenz 04 model
//' @export
// [[Rcpp::export]]
arma::vec RHS(const arma::vec& X, const int& K, const double& F){
  
  arma::vec W = getW(X, K);
  arma::vec XXnew = getXX(X, K, W);
  arma::vec newX = XXnew - X + F;
 
  return(newX); 
}
