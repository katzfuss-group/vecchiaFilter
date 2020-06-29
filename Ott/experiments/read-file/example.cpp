#include "vector.h"
#include <iostream>


int main(){

  for( const double& i : v)
    std::cout << i << ' ';
  std::cout << '\n';
  
  return 0;
  
}
