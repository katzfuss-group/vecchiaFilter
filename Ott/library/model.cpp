#include <cstdlib>
#include <cmath>
#include <fstream>
#include "math_lib.h"
#include "model.h"

Lorenz_model::Lorenz_model( int model_number, int N, double F, double dt, double b, double c )
: model_number( model_number ), N( N ), K( N/30 ), F( F ), dt( dt ), b( b ), c( c ), Z( N ), X( N ), Y( N ), WK( N ), bracket( N ), k1( N ), k2( N ), k3( N ), k4( N )
{
	I = set_I();
	alpha = double(3*I*I+3) / (2*I*I*I+4*I); 
	beta = double(2*I*I+1) / (I*I*I*I+2*I*I);
	randomize();
}

int Lorenz_model::size() const { return N; }

void Lorenz_model::randomize() { for( int i=0; i<N; ++i ) Z[i] = random_uniform(); }

void Lorenz_model::set_state( const vector &state )
{
	for( int d=state.size()/N, n=0; n<N; ++n ) Z[n] = state[n*d];
}

void Lorenz_model::set_state( double *state, int size )
{
	for( int d=size/N, n=0; n<N; ++n ) Z[n] = state[n*d];
}

void Lorenz_model::get_state( vector &state ) const
{
	for( int d=N/state.size(), n=0; n<state.size(); ++n ) state[n] = Z[n*d];
}

void Lorenz_model::get_state( double *state ) const
{
	for( int i=0; i<N; ++i ) state[i] = Z[i];
}

const vector &Lorenz_model::get_state()
{
	return Z;
}

double *Lorenz_model::get_v()
{
	return Z.get_v();
}

void Lorenz_model::get_XandY( vector &X, vector &Y )
{
	make_XandY( Z, X, Y ); 
}

void Lorenz_model::evolve( double time )
{
	// vector k1( N ), k2( N ), k3( N ), k4( N );

	int steps = int( time/dt + 0.5 );

	for( int i = 0; i < steps; ++i )
	{
		dZ( Z,         k1 );
		dZ( Z + k1/2., k2 );
		dZ( Z + k2/2., k3 );
		dZ( Z + k3,    k4 );
		Z += ( k1 + 2.*k2 + 2.*k3 + k4 ) / 6.;
	}
}

// private member functions

inline int Lorenz_model::mod( int i )
{
	int j = i%N;
	return j<0 ? j+N : j;
}

void Lorenz_model::dZ( const vector &Z, vector &dZ )
{
	switch( model_number )
	{
		case 1: { dZ_1( Z, dZ ); return; }
		case 2: { dZ_2( Z, dZ ); return; }
		case 3: { dZ_3( Z, dZ ); return; }
	}
}

void Lorenz_model::dZ_1( const vector &Z, vector &dZ )
{
	for( int n=0; n<N; ++n ) dZ[n] = dt * ( -Z[ mod(n-2) ] * Z[ mod(n-1) ] + Z[ mod(n-1) ] * Z[ mod(n+1) ] - Z[n] + F );
}

void Lorenz_model::dZ_2( const vector &Z, vector &dZ )
{
	// vector bracket( N );
	
	make_bracket( Z, bracket );
	for( int n=0; n<N; ++n ) dZ[n] = dt * ( bracket[n] - Z[n] + F );
}

void Lorenz_model::dZ_3( const vector &Z, vector &dZ )
{
	// vector X( N ), Y( N ), bracket( N );

	make_XandY( Z, X, Y );
	make_bracket( X, bracket );
	for( int n=0; n<N; ++n )
	{
		dZ[n] = dt * (
		b * b * ( -Y[ mod(n-2) ] * Y[ mod(n-1) ] + Y[ mod(n-1) ] * Y[ mod(n+1) ] )
		+ c * ( -Y[ mod(n-2) ] * X[ mod(n-1) ] + Y[ mod(n-1) ] * X[ mod(n+1) ] ) 
		+ bracket[n] - X[n] - b * Y[n] + F );
	}
}

void Lorenz_model::make_XandY( const vector &Z, vector &X, vector &Y )
{
	double sum;

	for( int n=0; n<N; ++n )
	{
		sum = alpha * Z[n];
		for( int i=1; i<I; ++i )
		{
			sum += ( alpha - beta*i ) * Z[ mod(n-i) ];
			sum += ( alpha - beta*i ) * Z[ mod(n+i) ];
        }
		sum += ( alpha - beta*I ) * Z[ mod(n-I) ] / 2.;
		sum += ( alpha - beta*I ) * Z[ mod(n+I) ] / 2.;
		X[n] = sum;
		Y[n] = Z[n] - X[n];
	}
}

void Lorenz_model::make_WK( const vector &X, vector &WK )
{
	bool flag_even = false;
	int J = K/2;
	double sum = 0.;
	if( K%2 == 0 ) flag_even = true;

	// Calculate WK[0]
	for( int j = -J; j <= J; ++j ) sum += X[ mod(0+j) ];
	WK[0] = sum;
	if( flag_even ) WK[0] -= ( X[ mod(0-J) ] + X[ mod(0+J) ] ) / 2.;

	// Calculate WK[n] for each n other than 0
	for( int n = 1; n < N; ++n )
	{
		sum += ( -X[ mod(n-J-1) ] + X[ mod(n+J) ] );
		WK[n] = sum;
		if( flag_even ) WK[n] -= ( X[ mod(n-J) ] + X[ mod(n+J) ] ) / 2.;
	}
}

void Lorenz_model::make_bracket( const vector &X, vector &bracket )
{
	bool flag_even = false;
	int J = K/2;
	double sum = 0.;
	// vector WK( N );
	if( K%2 == 0 ) flag_even = true;

	make_WK( X, WK );
	
	// Calculate bracket[0]
	for( int j = -J; j <= J; ++j ) sum += WK[ mod(0-K+j) ] * X[ mod(0+K+j) ];
	bracket[0] = - WK[ mod(0-2*K) ] * WK[ mod(0-K) ] + sum;
	if( flag_even ) bracket[0] -= ( WK[ mod(0-K-J) ] * X[ mod(0+K-J) ] + WK[ mod(0-K+J) ] * X[ mod(0+K+J) ] ) / 2.;

	// Calculate bracket[n] for the rest of n's
	for( int n = 1; n < N; ++n )
	{
		sum += - WK[ mod(n-K-J-1) ] * X[ mod(n+K-J-1) ] + WK[ mod(n-K+J) ] * X[ mod(n+K+J) ] ;
		bracket[n] = - WK[ mod(n-2*K) ] * WK[ mod(n-K) ] + sum;
		if( flag_even ) bracket[n] -= ( WK[ mod(n-K-J) ] * X[ mod(n+K-J) ] + WK[ mod(n-K+J) ] * X[ mod(n+K+J) ] ) / 2.;
	}
	
	bracket /= (K*K);
}

int Lorenz_model::set_I() // make 7 waves, maximum I values
{
	switch( N )
	{
		case 960: return 12;
		case 480: return 6;
		case 240: return 4;
		case 120: return 3;
		case  60: return 1;
		case  30: return 1;
	}
}
