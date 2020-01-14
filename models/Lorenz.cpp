#include <iostream>
#include <armadillo>



int mod_cpp(const int i, const int n){ 
  int mod = (i - 1) % n;
  if(mod < 0)  {
    mod += n;
  }
  return mod;
}

double Wn_Even_Cpp(const arma::colvec& XX, const int& n, const int& k, const int& N)
{
  double K = k;
  double J = K / 2;
  
  double first = XX(mod_cpp(n - J, N)) / (2 * K);
  double last = XX(mod_cpp(n + J, N)) / (2 * K);
  
  double sums = first + last - XX(mod_cpp(n, N)) / K;
  
  for(int i = 0; i < J; ++i)
  {
    sums = sums + XX(mod_cpp(n - i, N)) / K + XX(mod_cpp(n + i, N)) / K;
  }
  
  return(sums);
}


double Wn_Odd_Cpp(const arma::colvec& XX, const int& n, const int& k, const int& N){
  double K = k;
  double J = (K - 1) /2;
  double sums = - XX(mod_cpp(n, N)) / K;
  
  for(int i = 0; i <= J; ++i)
  {
    sums = sums + XX(mod_cpp(n - i, N)) / K + XX(mod_cpp(n + i, N)) / K;
  }
  
  return(sums);
}


double XX_Kn_Even_Cpp(const arma::colvec& XX, const int& n, const int& k, 
                      const int& N)
{
  double K = k;
  double J = K / 2;
  
  double first = (Wn_Even_Cpp(XX, mod_cpp(n - K - J, N), K, N) * XX(mod_cpp(n + K - J, N))) / (2 * K);
  double last = (Wn_Even_Cpp(XX, mod_cpp(n - K + J, N), K, N) * XX(mod_cpp(n + K + J, N))) / (2 * K);
  
  double r_sum = first + last + (Wn_Even_Cpp(XX, mod_cpp(n - K, N), K, N) * XX(mod_cpp(n + K, N)) / K);
  
  for(int i = 1; i < J; ++i){
    r_sum = r_sum + Wn_Even_Cpp(XX, mod_cpp(n - K - i, N), K, N) * XX(mod_cpp(n + K - i, N)) / K +
            Wn_Even_Cpp(XX, mod_cpp(n - K + i, N), K, N) * XX(mod_cpp(n + K + i, N)) / K;
    //std::cout << "iteration " << i << ": " << r_sum << std::endl;
  }
  
  double XX_Kn_val = (-Wn_Even_Cpp(XX, mod_cpp(n - (2 * K), N), K, N)) * Wn_Even_Cpp(XX, mod_cpp(n - K, N), K, N) + r_sum;
  //std::cout << (-Wn_Even_Cpp(XX, mod_cpp(n - (2 * K), N), K, N)) * Wn_Even_Cpp(XX, mod_cpp(n - K, N), K, N) << std::endl;
  return XX_Kn_val;
}


double XX_Kn_Odd_Cpp(const arma::colvec& XX, const int& n, const int& k, const int& N){
  
  double K = k;
  double J = (K - 1) / 2;
  
  double r_sum = (Wn_Odd_Cpp(XX, mod_cpp(n - K, N), K, N) * XX(mod_cpp(n + K, N)) / K);
  
  for(int i = 1; i <= J; ++i){
    r_sum = r_sum + Wn_Odd_Cpp(XX, mod_cpp(n - K - i, N), K, N) * XX(mod_cpp(n + K - i, N)) / K +
      Wn_Odd_Cpp(XX, mod_cpp(n - K + i, N), K, N) * XX(mod_cpp(n + K + i, N)) / K;
    //std::cout << "iteration " << i << ": " << r_sum << std::endl;
  }
  
  double XX_Kn_val = (-Wn_Odd_Cpp(XX, mod_cpp(n - (2 * K), N), K, N)) * Wn_Odd_Cpp(XX, mod_cpp(n - K, N), K, N) + r_sum;
  
  return XX_Kn_val;
  
}


arma::colvec Lorenz04M2Cpp(const arma::colvec& X0, const int & F_Lor, const int& K_Lor,
                           const double& dt, const int& M, const int& N_Lor)
{
  // Create vector to return with passed vector
  arma::colvec XX = X0;
  // Initialize vector of zeros that will "update" XX
  arma::colvec dx(N_Lor);
  dx.fill(0);
  double F = F_Lor;
  
  if(K_Lor % 2 == 0){
    // Nested for loops to update XX
    for(int i = 0; i < M ; ++i){
      for(int j = 0; j < N_Lor; ++j){
        dx(j) = (XX_Kn_Even_Cpp(XX, j, K_Lor, N_Lor) - XX(j) + F) * dt;
        //std::cout << dx << std::endl;
      }
      
      XX = XX + dx;
      //std::cout << XX << std::endl;
      
      
    }
  }
  else
  {
    // Nested for loops to update XX
    for(int i = 0; i < M ; ++i){
      for(int j = 0; j < N_Lor; ++j){
        dx(j) = (XX_Kn_Odd_Cpp(XX, j, K_Lor, N_Lor) - XX(j) + F) * dt;
        //std::cout << dx << std::endl;
      }
      
      XX = XX + dx;
      //std::cout << XX << std::endl;
      
      
    }
  }
  
  
  
  return XX;
}



arma::mat Lorenz04M2SimCpp(const arma::colvec& Xinit, const int& F_Lor, const int& K_Lor,
                        const double& dt, const int& M, const int& iter, const int& burn)
{
  // Get N_Lor from Xinit
  int N_Lor = Xinit.n_rows;
  // Create Matrix to store generated data
  arma::mat Xiter(N_Lor, iter);
  
  // Create vector to allow for burn-in iterations
  arma::colvec newXburn = Xinit;
  
  int burnin = burn;
  
  for(int i = 0; i < burnin; ++i){
    newXburn = Lorenz04M2Cpp(newXburn, F_Lor, K_Lor, dt, M, N_Lor);
  }
  
  Xiter.col(0) = Lorenz04M2Cpp(newXburn, F_Lor, K_Lor, dt, M, N_Lor);
  
  for(int i = 1; i < iter; ++i){
    Xiter.col(i) = Lorenz04M2Cpp(Xiter.col(i - 1), F_Lor, K_Lor, dt, M, N_Lor);
  }
  
  
  return Xiter;
}




int main(int argc, const char **argv) {

    arma::arma_rng::set_seed_random();

    int N = 4;
    int F = 10;
    int K = 32;
    double dt = 0.005;
    int M = 40;
    arma::colvec X0 = arma::randn(N);

    arma::mat X = Lorenz04M2SimCpp(X0, F, K, dt, M, 2,0);

    std::cout << X << std::endl;
      
    return 0;
}

