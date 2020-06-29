#include <algorithm>
#include <iostream>
#include <fstream>
#include "math_lib.h"
#include "observation.h"

observation::observation( double sigma, const std::vector<int> &locations, int size_true_state ) 
: sigma(sigma), n_ob(locations.size()), N(size_true_state), location(0), mark(0)
{
	set_locations( locations );	
}

observation::~observation() { delete [] location; delete [] mark; }

void observation::set_locations( const std::vector<int> &locations )
{
	delete [] location;
	delete [] mark;

	n_ob = locations.size();
	location = new int[n_ob];
	ob.resize(n_ob);
	mark = new int[N];

	for( int i=0; i<N; ++i ) mark[i] = -1;
	for( int i=0; i<n_ob; ++i ) location[i] = locations[i];	
	for(int i=0; i<n_ob; ++i) { mark[location[i]] = i; }
}

void observation::state_to_ob(const vector &ApprState, vector &ob) const
{
	int Na = ApprState.size();
	int d = N/Na;
	ob.resize(n_ob);

	for(int na0, na1, i=0; i<n_ob; ++i)
	{
		na0 = location[i]/d;
		na1 = ( na0+1<Na ? na0+1 : 0 );

		ob[i] = ApprState[na0]
			 + static_cast<double>( (ApprState[na1] - ApprState[na0] ) * ( location[i]-na0*d ) ) / d;
	}
}

void observation::Jacobian(const vector &ApprState, matrix &H) const
// calculate global Jacobian matrix of observation operator
{
	int Na = ApprState.size();
	int d = N/Na;
	H.resize(n_ob, Na);
	for(int i=0; i<n_ob; ++i) for(int j=0; j<Na; ++j) H[i][j] = 0;
	
	for(int na0, na1, i=0; i<n_ob; ++i)
	{
		na0 = location[i]/d;
		na1 = ( na0+1<Na ? na0+1 : 0 );
		
		H[i][na0] = 1. - double( location[i] - na0*d ) / d;
		H[i][na1] = 1. - H[i][na0];
	}
}

/*void observation::state_to_ob(const vector &ApprState, vector &ob) const
{
	int Na = ApprState.size();
	int d = N/ApprState.size();
	int d3 = d*d*d;
	ob.resize(n_ob);
	
	int n0, n1, n2, n3;
	double z0, z1, z2, z3;
	int Lag0, Lag1, Lag2, Lag3;
	
	Lag0 = -6*d3; Lag1 = 2*d3; Lag2 = -2*d3; Lag3 = 6*d3;
	
	for(int i=0; i<Na; ++i)
	{
		n0 = i*d;
		n1 = n0 + d;
		n2 = n1 + d;
		n3 = n2 + d;
		
		z0 = ApprState[i];
		z1 = ApprState[(i+1)%Na];
		z2 = ApprState[(i+2)%Na];
		z3 = ApprState[(i+3)%Na];
		
		//Lag0 = (n0-n1)*(n0-n2)*(n0-n3);
		//Lag1 = (n1-n0)*(n1-n2)*(n1-n3);
		//Lag2 = (n2-n0)*(n2-n1)*(n2-n3);
		//Lag3 = (n3-n0)*(n3-n1)*(n3-n2);
		
		if(mark[n0] > -1) ob[mark[n0]] = z0;
		
		for(int j, n=n1+1; n<n2; ++n)
		{
			j = mark[n%N];
			if (j > -1) ob[j] = z0*(n-n1)*(n-n2)*(n-n3)/Lag0
					+ z1*(n-n0)*(n-n2)*(n-n3)/Lag1
					+ z2*(n-n0)*(n-n1)*(n-n3)/Lag2
					+ z3*(n-n0)*(n-n1)*(n-n2)/Lag3;
		}
	}
}*/

void observation::ob_to_state(vector &ApprState) const
{
	int d = N/ApprState.size();
	int n0, n1;
	double z0, z1;

	for(int i=0; i<n_ob; ++i)
	{
		n0 = location[i];
		n1 = location[(i+1)%n_ob];
		z0 = ob[i];
		z1 = ob[(i+1)%n_ob];
		
		n1 = n0 + cycle(n1-n0, N);
		
		for( int n=(n0+d-1)/d*d; n<n1; n+=d) {
			ApprState[cycle(n,N)/d] = z0 + double( (z1-z0)*(n-n0) ) / (n1-n0);
		}
	}
}

/*void observation::ob_to_state(vector &ApprState) const
{
	int d = N/ApprState.size();
	int n0, n1, n2, n3;
	double z0, z1, z2, z3;
	int Lag0, Lag1, Lag2, Lag3;

	for(int i=0; i<n_ob; ++i)
	{
		n0 = location[i];
		n1 = location[(i+1)%n_ob];
		n2 = location[(i+2)%n_ob];
		n3 = location[(i+3)%n_ob];

		z0 = ob[i];
		z1 = ob[(i+1)%n_ob];
		z2 = ob[(i+2)%n_ob];
		z3 = ob[(i+3)%n_ob];

		n1 = n0 + cycle(n1-n0, N);
		n2 = n0 + cycle(n2-n0, N);
		n3 = n0 + cycle(n3-n0, N);
	
		Lag0 = (n0-n1)*(n0-n2)*(n0-n3);
		Lag1 = (n1-n0)*(n1-n2)*(n1-n3);
		Lag2 = (n2-n0)*(n2-n1)*(n2-n3);
		Lag3 = (n3-n0)*(n3-n1)*(n3-n2);

		for( int n=(n1+d-1)/d*d; n<n2; n+=d) {
			ApprState[cycle(n,N)/d] = static_cast<double>(z0*(n-n1)*(n-n2)*(n-n3))/Lag0
						+ static_cast<double>(z1*(n-n0)*(n-n2)*(n-n3))/Lag1
						+ static_cast<double>(z2*(n-n0)*(n-n1)*(n-n3))/Lag2
						+ static_cast<double>(z3*(n-n0)*(n-n1)*(n-n2))/Lag3;
		}
	}
}*/

void observation::cov_inv( int middle, int span, int size_ApprState, vector &r_inv ) const
{
	int d = N / size_ApprState;
	int middle_true = d * middle; // location of middle in the true state
	int span_true = d * span; // span in the true state
	int r = span_true * 5 / 5, s = span_true - r, n_ob_l = 0; // number of local observations in the local region

	for( int l = -span_true; l <= span_true; ++l ) if( mark[ cycle(middle_true+l,N) ] > -1 ) ++n_ob_l;

	r_inv.resize( n_ob_l );

	for( int i=0, l = -span_true; l <= span_true; ++l ) if( mark[ cycle(middle_true+l,N) ] > -1 )
	{
		if( abs(l) > r ) r_inv[i++] = ( 1. / ( sigma * sigma ) ) * ( 1. - double(abs(l)-r)/s );
		else r_inv[i++] = 1. / ( sigma * sigma );
	}
}

void observation::global_to_local( int middle, int span, int size_ApprState, const vector &ob_g, vector &ob_l ) const
{
	// distinguish between a location and an index of location
	// local region refers to a local region in the true state
	int d = N / size_ApprState;
	int middle_true = d * middle; // location of middle in the true state
	int span_true = d * span; // span in the true state
	int l1 = cycle(middle_true - span_true, N); // leftmost location of local region
	int n_ob_l = 0; // number of local observations in the local region

	for(int l=0; l<2*span_true+1; ++l) if( mark[cycle(l1+l, N)] > -1 ) ++n_ob_l;

	ob_l.resize(n_ob_l);

	for(int m, i=0, l=0; l<2*span_true+1; ++l) {
		m = mark[cycle(l1+l, N)];	
		if( m > -1 ) ob_l[i++] = ob_g[m];
	}
}

void observation::global_to_local_Jacobian
(int middle, int span, int size_ApprState, const matrix &Hg, matrix &H) const
{
	int d = N / size_ApprState;
	int middle_true = d * middle;
	int span_true = d * span;
	int l1 = cycle(middle_true - span_true, N);
	int j1 = cycle(middle - span, size_ApprState);
	int n_ob_l = 0;

	for(int l=0; l<2*span_true+1; ++l) if( mark[cycle(l1+l, N)] > -1 ) ++n_ob_l;

	H.resize(n_ob_l, 2*span+1);

	for(int m, i=0, l=0; l<2*span_true+1; ++l) {
		m = mark[cycle(l1+l, N)];
		if( m > -1 ) {
			for(int j=0; j<2*span+1; ++j) {
				H[i][j] = Hg[m][cycle(j1+j, size_ApprState)];
			}
			++i;
		}
	}
}

void observation::global_to_region( int region_begin, int region_end, const vector &ob_g, vector &ob_r ) const
{
	int n_ob_r = 0; // number of local observations in the local region

	for( int n = region_begin; n <= region_end; ++n ) if( mark[ cycle( n,N ) ] > -1 ) ++n_ob_r;

	ob_r.resize( n_ob_r );

	for( int m, i = 0, n = region_begin; n <= region_end; ++n )
	{
		m = mark[ cycle( n, N ) ];
		if( m > -1 ) ob_r[i++] = ob_g[m];
	}
}

std::ostream& operator << (std::ostream& os, const observation &obsv)
{
	vector ob;
	ob = obsv.get_observation();	
	os << ob;
	return os;
}

