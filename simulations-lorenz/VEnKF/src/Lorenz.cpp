#include "oldRHS.h"
#include "newRHStools.h"
#include <RcppArmadillo.h>
#include "exactGrad.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;





//' A single iteration of the Lorenz 04 model
//'
//' @param XX initial state
//' @param F_Lor F from the Lorenz 04 model
//' @param K_Lor K from the Lorenz 04 model
//' @param dt time step
//' @param M the number of steps that make a time step
//' @param order the order of RK method to use
//' @return dx the increment to reach the new state
//' @export
// [[Rcpp::export]]
arma::vec DeltaLorenz04M2Cpp(const arma::vec X0, const double & F, const int& K,
                             const double& dt, const int& M, const bool& newAlgo){
  
  // Initialize vector of zeros that will "update" XX
  arma::vec X = X0;
  int N_Lor = X.n_rows;

  arma::vec dx;
  
  // Nested for loops to update XX
  for(int i = 0; i < M ; ++i){

    if(newAlgo){
      
      arma::vec k1 = RHS(X, K, F);
      arma::vec k2 = RHS(X + 0.5*dt*k1, K, F);
      arma::vec k3 = RHS(X + 0.5*dt*k2, K, F);
      arma::vec k4 = RHS(X + dt*k3, K, F);
      
      dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
      
    } else {

      dx = arma::zeros<arma::vec>(N_Lor);
      for(int j = 0; j < N_Lor; ++j){
      
	      double k1 = oldRHS(X, j, K, F);
	      double k2 = oldRHS(X + 0.5*dt*k1, j, K, F);
	      double k3 = oldRHS(X + 0.5*dt*k2, j, K, F);
	      double k4 = oldRHS(X + dt*k3, j, K, F);
	      dx(j) =  (k1 + 2*k2 + 2*k3 + k4) * dt/6;
  
      }
    } 
    
    X += dx;
      
  }
  
  return X - X0;
}


//' Simulate from the Lorenz 04 model
//'
//' @param Xinit initial state
//' @param F_Lor F from the Lorenz 04 model
//' @param K_Lor K from the Lorenz 04 model
//' @param dt time step
//' @param M the number of steps that make a time step
//' @param iter number of iterations
//' @param burn how many steps to discard
//' @return Xiter values of all iterations
//' @export
// [[Rcpp::export]]
arma::mat Lorenz04M2SimCpp(const arma::vec& Xinit, const int& F_Lor, const int& K_Lor,
                           const double& dt, const int& M, const int& iter, const int& burn,
			   const bool& newAlgo)
{
  // Get N_Lor from Xinit
  int N_Lor = Xinit.n_rows;
  // Create Matrix to store generated data
  arma::mat Xiter(N_Lor, iter+1);
  
  // Create vector to allow for burn-in iterations
  arma::vec newXburn = Xinit;
  
  int burnin = burn;
  
  if( burnin > 0 ){
    for(int i = 0; i < burnin; ++i){
      arma::vec delta = DeltaLorenz04M2Cpp(newXburn, F_Lor, K_Lor, dt, M, newAlgo);
      newXburn = newXburn + delta;
    }
  } 
  
  Xiter.col(0) = newXburn;
  
  for(int i = 1; i < iter+1; ++i){
    arma::vec delta = DeltaLorenz04M2Cpp(Xiter.col(i - 1), F_Lor, K_Lor, dt, M, newAlgo);
    Xiter.col(i) = Xiter.col(i-1) + delta;
  }
  
  return Xiter.tail_cols(iter);
}





arma::mat EGradient(const arma::vec& X, const int& K, const double& dt, const double& F){
  
  int N = X.n_rows;
  
  arma::vec Xp = X + 0.5*dt*RHS(X, K, F);
  arma::vec Xb = X + 0.5*dt*RHS(Xp, K, F);
  arma::vec Xt = X + dt*RHS(Xb, K, F);
  
  arma::mat A = exactGrad(X, K);
  arma::mat B = exactGrad(Xp, K) % (1 + 0.5*dt*A);
  arma::mat C = exactGrad(Xb, K) % (1 + 0.5*dt*B);
  arma::mat D = exactGrad(Xt, K) % (1 + dt*C);
  
  return (A + B + C + D)*(dt/6) + arma::eye(N,N);
}



//' Gradient for the Lorenz model
//'
//' @param Xinit initial state
//' @param K_Lor K from the Lorenz 04 model
//' @param M the number of steps that make a time step
//' @param dt time step
//' @param F_Lor F from the Lorenz 04 model
//' @export
// [[Rcpp::export]]
arma::mat exactGradient(const arma::vec& X, const int& K, const int& M, const double& dt, const double& F){

  int N = X.n_rows;
  arma::mat G = arma::eye(N, N);
  arma::vec Xnew = X;
  for( int m=1; m<=M; m++){
    G = EGradient(Xnew, K, dt, F) * G;
    Xnew = Xnew + DeltaLorenz04M2Cpp(Xnew, F, K, dt, 1, false);
  }
  return( G );
}




arma::mat Lorenz04M2Gradient(const arma::vec& X, const int& F_Lor, const int& K_Lor,
                           const double& dt, const int& M, const int& iter, const int& burn,
			     const bool& newAlgo){

  
  int N = X.n_rows;
  arma::vec newf(N);
  arma::vec del(N);
  arma::vec tempX = X;
  arma::vec reff = X + DeltaLorenz04M2Cpp(X, F_Lor, K_Lor, dt, M, newAlgo);
  arma::mat gradient = arma::zeros<arma::mat>(N,N);

  arma::vec pert = arma::ones<arma::vec>(N)*1e-8;
  arma::vec delta = arma::max(pert, arma::abs(X % pert));
  
  
  for(arma::uword i=0; i<N; i++){
  
    tempX(i) = tempX(i) + delta(i);
    newf = tempX + DeltaLorenz04M2Cpp(tempX, F_Lor, K_Lor, dt, M, newAlgo);
    del = (newf - reff) / delta(i);
    
    gradient.col(i) = del;
    tempX(i) = X(i);
  }
  
  return gradient;

}
