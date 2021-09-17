#include <RcppArmadillo.h>
#include "vectorTools.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



double sumPrime(arma::vec v, const int& K){
  if( K % 2==0 ){
    v.head(1) /=2;
    if(v.size()>1){
      v.tail(1) /= 2;
    }
  }
  return(sum(v));
}


bool check(int ind1, int ind2, int p){
  if(ind1 < ind2){
    return p>=ind1 && p<=ind2;
  } else if(ind1==ind2) {
    return p==ind1;
  } else {
    return p>=ind1 || p<=ind2;
  }
}


bool isIn1( int n, int p, int K, int J, int N ){
  int ind1 = mod(n - 2*K - J, N);
  int ind2 = mod(n - 2*K + J, N);
  return check(ind1, ind2, p);
}


bool isIn2( int n, int p, int K, int J, int N ){
  int ind1 = mod(n - K - J, N);
  int ind2 = mod(n - K + J, N);
  
  return check(ind1, ind2, p);
}



bool isIn3( int n, int p, int K, int J, int N ){
  int ind1 = mod(n - K - 2*J, N);
  int ind2 = mod(n - K + 2*J, N);
  // if(p==19 && n==0){
  //   Rcout << "n=" << n << std::endl;
  //   Rcout << "p=" << p << std::endl;
  //   Rcout << "K=" << K << std::endl;
  //   Rcout << "J=" << J << std::endl;
  //   Rcout << "N=" << N << std::endl;
  //   Rcout << "ind1=" << ind1 << std::endl;
  //   Rcout << "ind2=" << ind2 << std::endl;
  // }
  return check(ind1, ind2, p);
}


bool isIn4( int n, int p, int K, int J, int N ){
  int ind1 = mod(n + K - J, N);
  int ind2 = mod(n + K + J, N);
  return check(ind1, ind2, p);
}



//' Calculate the exact gradient in the Lorenz model
//'
//' @param X the vector at which to evaluate the gradient
//' @param K parameter from the Lorenz model
//' @return grad the matrix with the gradient
//' @export
// [[Rcpp::export]]
arma::mat exactGrad(const arma::vec& X, const int& K) {
  
  int N = X.n_rows;
  mat grad(N,N);
  int J = K/2;
  //Rcout << "J=" << J << std::endl;
  
  uword ind1, ind2;
  
  for( int n = 0; n<N; n++ ){
    //Rcout << "======= n=" << n << " =======" << std::endl;
    for( int p = 0; p<N; p++ ){
      //Rcout << "=== p=" << p << "===" << std::endl;
      double d = 0;
      if( isIn1(n, p, K, J, N) ){
        //Rcout << "in 1" << std::endl;
        ind1 = mod(n - K - J, N);
        ind2 = mod(n - K + J, N);
        vec masked = mask(X, ind1, ind2, true);
        d -= sumPrime(masked, K);
      }
      if( isIn2(n, p, K, J, N) ){
        //Rcout << "in 2" << std::endl;
        ind1 = mod(n - 2*K - J, N);
        ind2 = mod(n - 2*K + J, N);
        vec masked = mask(X, ind1, ind2, true);
        d -= sumPrime(masked, K);
      }
      if( isIn3(n, p, K, J, N) ){
        //Rcout << "in 3" << std::endl;
        for( int j=-J; j<=J; j++){
          for( int i=-J; i<=J; i++){
            if(mod(n - K + j - i, N) == p){
              d += X(mod(n + K + j, N)) * (abs(j)==J)/2;      
            }
          }
        }
      }
      if( isIn4(n, p, K, J, N) ){
        //Rcout << "in 4" << std::endl;
        int j = p - n - K;
        ind1 = mod(n - K + j - J, N);
        ind2 = mod(n - K + j + J, N);
        vec masked = mask(X, ind1, ind2, true);
        d -= sumPrime(masked, K);
      }
      if( p==n ){
        d -= 1;
      }
      grad(n, p) = d/(K*K);
    }
  }
  return grad;
}

