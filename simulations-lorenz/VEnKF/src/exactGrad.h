#include <RcppArmadillo.h>
#include "vectorTools.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



bool check(int ind1, int ind2, int p);

bool isIn1( int n, int p, int K, int J, int N );

bool isIn2( int n, int p, int K, int J, int N );

bool isIn3( int n, int p, int K, int J, int N );

bool isIn4( int n, int p, int K, int J, int N );

mat exactGrad(const vec& X, const int& K);