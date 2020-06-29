#include <iostream>
#include "tensor.h"
#include "model.h"
#include "initval.h"

using namespace std;

int main(){

  int model_number = 3;
  int N = 960;
  double F = 15;
  double dt = 0.005;
  double b = 10;
  double c = 2.5;
  double Tmax = 100*dt;
  
  vector x0vec(N);
  
  for( int i = 0; i < N; ++i ){
    x0vec[i] = x0[i];
  }

  LM L = LM( model_number, N, F, dt, b, c);

  cout << "Z_0: ";
  for( int i = 0; i<3; ++i ){
    cout << x0vec[i] << ' ';
  }
  cout << endl;
  
  
  L.set_state( x0vec );
  L.evolve(Tmax);
  
  vector Z = L.get_state();

  cout << "Z_1: ";
  for( int i = 0; i<3; ++i ){
    cout << Z[i] << ' ';
  }
  cout << endl;
  
}

