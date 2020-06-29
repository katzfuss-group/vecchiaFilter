#ifndef MATH_LIB_H
#define MATH_LIB_H

#define LAPACK
#include <vector>
#include "tensor.h"

// random number generators
double random_uniform();
double random_gaussian( double sigma );

// matrix functions
void mul( const matrix &a, const matrix &b, matrix &c );
matrix operator * ( const matrix &a, const matrix &b );
int inverse( matrix &a );
int inverse( const matrix &a, matrix &aInv );
int eigen( matrix &a, vector &eigenvalues );
int eigen( const matrix &a, vector &eigenvalues, matrix &eigenvectors );

// other
void DFT( const vector &x, vector &real, vector &imag, int kmin = 1, int kmax = 0 );
void DFT( const matrix &z, matrix &real, matrix &imag, int kxmin = 1, int kxmax = 0, int kymin = 1, int kymax = 0 );
void envelope( const vector &RealFT, const vector &ImagFT, vector &env, int kmin = 1, int kmax = 0 );
void envelope( const matrix &RealFT, const matrix &ImagFT, matrix &env, 
	int kxmin = 1, int kxmax = 0, int kymin = 1, int kymax = 0 );
void choose_regular_numbers( int n, int N, std::vector<int> &numbers ); // choose n numbers from 0 to N-1
void choose_random_numbers( int n, int N, std::vector<int> &numbers ); // choose n numbers from 0 to N-1 
int downhill_simplex( double (*func)( const vector & ), matrix &vertex, double termination );
double find_root( double (*func)( double ), double (*funcd)( double ), double x1, double x2, double accuracy );
inline int mod( const int &n, const int &N ) { int m = n % N; return m < 0 ? m + N : m; }

// numerical recipe
namespace NR
{
	template<class T> inline const T SQR(const T a) {return a*a; }
	template<class T> inline const T SIGN(const T &a, const T &b)
		{ return b>=0 ? ( a>=0 ? a : -a ) : ( a>=0 ? -a : a ); }
	template<class T> inline void SWAP(T &a, T &b) { T dum = a; a = b; b = dum; }
	template<class T> inline const T MAX( const T &a, const T &b ) { return b > a ? (b) : (a); }

	double pythag(const double &a, const double &b);
	int gaussj(matrix &a, matrix &b);
	void tred2(matrix &a, vector &d, vector &e);
	void tred2_trans(matrix &a, vector &d, vector &e);
	int tqli(vector &d, vector &e, matrix &z);
	int tqli_trans(vector &d, vector &e, matrix &z);
	double ran0( int &idum );
	double ran1( int &idum );
	void mnbrak( double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, double func( const double ) );
	double brent
	( const double ax, const double bx, const double cx, double f( const double ), const double tol, double &xmin );
	void linmin( vector &p, vector &xi, double &fret, double func( const vector & ) );
	double f1dim( const double x );
	void frprmn( vector &p, const double ftol, int &iter, double &fret, double func( const vector & ), 
		void dfunc( const vector &, vector & ) );
}

#endif
