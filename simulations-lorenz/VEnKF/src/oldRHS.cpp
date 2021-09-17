#include "vectorTools.h"
#include "newRHStools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double Wn_Even_Cpp(const arma::colvec& XX, const int& n, const int& K, const int& N){
  
  double J = K / 2;

  double first = XX(mod(n - J, N)) / (2 * K);
  double last = XX(mod(n + J, N)) / (2 * K);
  double sums = first + last - XX(mod(n, N)) / K;
  
  for(int i = 0; i < J; ++i)
  {
    sums = sums + XX(mod(n - i, N)) / K + XX(mod(n + i, N)) / K;
  }

  return(sums);
}


double Wn_Odd_Cpp(const arma::colvec& XX, const int& n, const int& k, const int& N){
  double K = k;
  double J = (K - 1) /2;

  double sums = - XX(mod(n, N)) / K;
  
  for(int i = 0; i <= J; ++i)
  {
    sums = sums + XX(mod(n - i, N)) / K + XX(mod(n + i, N)) / K;
  }

  return(sums);
}


double XX_Kn_Even_Cpp(const arma::colvec& XX, const int& n, const int& k){
  
  int N = XX.n_rows;
  double K = k;
  double J = K / 2;
  
  double first = (Wn_Even_Cpp(XX, mod(n - K - J, N), K, N) * XX(mod(n + K - J, N))) / (2 * K);
  double last = (Wn_Even_Cpp(XX, mod(n - K + J, N), K, N) * XX(mod(n + K + J, N))) / (2 * K);
  
  double r_sum = first + last + (Wn_Even_Cpp(XX, mod(n - K, N), K, N) * XX(mod(n + K, N)) / K);
  
  for(int i = 1; i < J; ++i){
    r_sum = r_sum + Wn_Even_Cpp(XX, mod(n - K - i, N), K, N) * XX(mod(n + K - i, N)) / K + Wn_Even_Cpp(XX, mod(n - K + i, N), K, N) * XX(mod(n + K + i, N)) / K;
    //std::cout << "iteration " << i << ": " << r_sum << std::endl;
  }
  
  double XX_Kn_val = (-Wn_Even_Cpp(XX, mod(n - (2 * K), N), K, N)) * Wn_Even_Cpp(XX, mod(n - K, N), K, N) + r_sum;
  //std::cout << (-Wn_Even_Cpp(XX, mod(n - (2 * K), N), K, N)) * Wn_Even_Cpp(XX, mod(n - K, N), K, N) << std::endl;
  return XX_Kn_val;
}


double XX_Kn_Odd_Cpp(const arma::colvec& XX, const int& n, const int& K){
  
  int N = XX.n_rows;
  //double K = k;
  double J = (K - 1) / 2;
  
  double r_sum = (Wn_Odd_Cpp(XX, mod(n - K, N), K, N) * XX(mod(n + K, N)) / K);
  
  for(int i = 1; i <= J; ++i){
    r_sum = r_sum + Wn_Odd_Cpp(XX, mod(n - K - i, N), K, N) * XX(mod(n + K - i, N)) / K + Wn_Odd_Cpp(XX, mod(n - K + i, N), K, N) * XX(mod(n + K + i, N)) / K;
    //std::cout << "iteration " << i << ": " << r_sum << std::endl;
  }
  
  double XX_Kn_val = (-Wn_Odd_Cpp(XX, mod(n - (2 * K), N), K, N)) * Wn_Odd_Cpp(XX, mod(n - K, N), K, N) + r_sum;
  
  return XX_Kn_val;
  
}


double oldRHS(const arma::vec& X, const int& j, const int& K, const double& F){

  // arma::vec W = getW(X, K);
  // Rcout << "hey" << std::endl;
  // arma::vec XXnew = getXX(X, K, W);
  // arma::vec newX = XXnew - X + F;
  // if( j==0 ){
  //   Rcout << "W(j)=" << W(j) << std::endl;
  //   Rcout << "Wn_Even_Cpp(j)=" << Wn_Even_Cpp(X, j, K, X.n_rows) << std::endl;
  // }
  
  double oldXj;
  double XXold;
  if(K % 2 == 0){
    XXold = XX_Kn_Even_Cpp(X, j, K);
    oldXj = XXold - X(j) + F;
  } else {
    oldXj = XX_Kn_Odd_Cpp(X, j, K) - X(j) + F;
  } 
  
  // Rcout << "diff=" << XXold - XXnew(j) << std::endl;
  // Rcout << "XXold=" << XXold << ", XXnew=" << XXnew(j) << std::endl;
  return(oldXj);
}
