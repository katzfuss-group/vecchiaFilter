#ifndef MODEL_H
#define MODEL_H

#include "tensor.h"

class Lorenz_model
{
	int model_number, N, K, I;
	double F, b, c, alpha, beta, dt; // dt = 0.05/7 is the maximum not to make model 3 blow up.
	vector Z, X, Y, WK, bracket, k1, k2, k3, k4;
	int mod( int i );
	int set_I();
	void make_XandY( const vector &Z, vector &X, vector &Y );
	void make_WK( const vector &X, vector &WK );
	void make_bracket( const vector &X, vector &bracket );
	void dZ( const vector &Z, vector &dZ );
	void dZ_1( const vector &Z, vector &dZ );
	void dZ_2( const vector &Z, vector &dZ );
	void dZ_3( const vector &Z, vector &dZ );
public:
	Lorenz_model( int model_number = 3, int N = 960, double F = 15., double dt = 0.05/12., double b = 10., double c = 2.5 );
	int size() const;
	void randomize();
	void set_state( const vector &state );
	void set_state( double *state, int size );
	void get_state( vector &state ) const;
	void get_state( double *state ) const;
	const vector &get_state();
	double *get_v();
	void get_XandY( vector &X, vector &Y );
	double &operator[] ( const int &index );
	const double &operator[] ( const int &index ) const;
	void evolve( double time );
};

inline double &Lorenz_model::operator[] (const int &index) { return Z[ index ]; }
inline const double &Lorenz_model::operator[] (const int &index) const { return Z[ index ]; }
typedef Lorenz_model LM;

#endif
