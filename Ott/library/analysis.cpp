#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "math_lib.h"
#include "analysis.h"

// Ensemble transform Kalman filter

ETKF::ETKF( double epsilon, char inflation_method ) : inflation_method( inflation_method )
{
	if( inflation_method == 'm' || inflation_method == 'x' ) inflation = 1. + epsilon;
	else if( inflation_method == 'r' ) inflation = epsilon;
	else { std::cerr << "not a valid inflation method\n"; exit( -1 ); }
}

int ETKF::analysis( matrix &X, const observation &obsv )
{
	int flag;
	//dimensions
	int k = X.size_row(); //number of ensemble members
	int m = X.size_col(); //dimension of the state
	int l = obsv.size(); //dimension of observation
	
	//for ensemble computation
	matrix Xa( k, m, 0 ); //analysis ensemble
	matrix Y( k); //observation ensemble corresponding to the state ensemble
	vector x_mean( m, 0 ); //mean of state ensemble
	vector y_mean( l, 0 ); //mean of observation ensemble

	// multiplicative inflation in state space
	if( inflation_method == 'x' )
	{
		vector x_mean_temp( m, 0. );
		for( int i = 0; i < k; ++i ) x_mean_temp += X[i];
		x_mean_temp /= k;
		
		for( int i = 0; i < k; ++i ) X[i] = inflation * ( X[i]-x_mean_temp ) + x_mean_temp;
	}
	
	//make Y and compute x_mean and y_mean
	for( int i=0; i<k; ++i )
	{
		obsv.state_to_ob( X[i], Y[i] );
		x_mean += X[i];
		y_mean += Y[i];
	}	
	x_mean /= k;
	y_mean /= k;

	//subtract x_mean and y_mean from X and Y
	for( int i=0; i<k; ++i )
	{
		X[i] -= x_mean;
		Y[i] -= y_mean;
	}
	
	vector yo = obsv.get_observation();
	matrix C;
	matrix Pw( k, k );
	matrix M( k, k ); //for eigenvectors
	matrix Minv( k, k );
	vector e( k ); //for eigenvalues
	matrix W( k, k );
	vector w_mean( k );
	vector temp( k );
	
	//compute C
	div( Y, obsv.covariance(), C );
	transpose( C );

	//compute Pw and W
	mul( Y, C, Pw );
	if( inflation_method == 'm' ) for( int i=0; i<k; ++i ) Pw[i][i] += (k-1) / inflation;
	else for( int i=0; i<k; ++i ) Pw[i][i] += (k-1);
	flag = eigen( Pw, e, M );
	if( flag == -1 ) return -1;
	inverse( M, Minv );
	if( flag == -2 ) return -2;
	for( int i=0; i<k; ++i ) for( int j=0; j<k; ++j ) M[i][j] *= sqrt( (k-1)/e[j] );
	mul( M, Minv, W );
	inverse( Pw );
	if( flag == -2 ) return -2;
	
	//compute w_mean
	yo -= y_mean;
	mul( yo, C, temp );
	mul( temp, Pw, w_mean );

	//relaxation
	if( inflation_method == 'r' )
	{
		mul( 1. - inflation, W );
		for( int i = 0; i < k; ++i ) W[i][i] += inflation;
	}
	
	//add w_mean to W
	for( int i=0; i<k; ++i ) W[i] += w_mean;

	//compute analysis ensemble
	mul( W, X, Xa );
	for( int i=0; i<k; ++i ) Xa[i] += x_mean;

	X = Xa;
	return 0;
}

// Ensemble transform Kalman filter for a sub-region

ETKF_regional::ETKF_regional( double epsilon, char inflation_method, int region_begin, int region_end )
: inflation_method( inflation_method ), region_begin( region_begin ), region_end( region_end )
{
	if( inflation_method == 'm' || inflation_method == 'x' ) inflation = 1. + epsilon;
	else if( inflation_method == 'r' ) inflation = epsilon;
	else { std::cerr << "not a valid inflation method\n"; exit( -1 ); }
}

int ETKF_regional::analysis( matrix &Xg, const observation &obsv )
{
	int flag;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = region_end - region_begin + 1; //dimension of local state
	
	matrix Yg( k ); //global observation ensemble corresponding to global state ensemble
	vector xg_mean( mg, 0 ); //mean of global state ensemble
	vector yg_mean( lg, 0 ); //mean of global observation ensemble

	// multiplicative inflation in state space
	if( inflation_method == 'x' )
	{
		vector x_mean_temp( mg, 0. );
		for( int i = 0; i < k; ++i ) x_mean_temp += Xg[i];
		x_mean_temp /= k;
		
		for( int i = 0; i < k; ++i ) Xg[i] = inflation * ( Xg[i]-x_mean_temp ) + x_mean_temp;
	}

	//make Yg and compute xg_mean and yg_mean
	for( int i = 0; i < k; ++i )
	{
		obsv.state_to_ob( Xg[i], Yg[i] );
		xg_mean += Xg[i];
		yg_mean += Yg[i];
	}	
	xg_mean /= k;
	yg_mean /= k;

	//subtract xg_mean and yg_mean from Xg and Yg
	for( int i = 0; i < k; ++i )
	{
		Xg[i] -= xg_mean;
		Yg[i] -= yg_mean;
	}
	
	//for regional computation
	matrix X( k, m ); //regional background ensemble
	matrix Y( k ); //regional observation ensemble
	vector x_mean( m );
	vector y_mean;
	vector yo; //regional observation
	matrix C;
	matrix Pw( k, k );
	matrix M( k, k ); //for eigenvectors
	matrix Minv( k, k );
	vector e( k ); //for eigenvalues
	matrix W( k, k );
	vector w_mean( k );
	vector temp( k );
	matrix Xa( k, m ); //regional analysis ensemble

	//get regional state and observation ensemble
	for( int i = 0; i < k; ++i )
	{
		get_region( Xg[i], X[i] );
		obsv.global_to_region( region_begin, region_end, Yg[i], Y[i]);
	}

	//get regional x_mean, y_mean, yo
	get_region( xg_mean, x_mean );
	obsv.global_to_region( region_begin, region_end, yg_mean, y_mean );
	obsv.global_to_region( region_begin, region_end, obsv.get_observation(), yo );

	if( yo.size() == 0 ) { Xa = X; goto final_Xa; }

	//compute C
	div( Y, obsv.covariance(), C );
	transpose(C);

	//compute Pw and W
	mul( Y, C, Pw );
	if( inflation_method == 'm' ) for( int i=0; i<k; ++i ) Pw[i][i] += (k-1) / inflation;
	else for( int i=0; i<k; ++i ) Pw[i][i] += (k-1);
	flag = eigen( Pw, e, M );
	if( flag == -1 ) return -1;
	flag = inverse( M, Minv );
	if( flag == -2 ) return -2;
	for( int i = 0; i < k; ++i ) for( int j=0; j < k; ++j ) M[i][j] *= sqrt( (k-1)/e[j] );
	mul( M, Minv, W );
	flag = inverse( Pw );
	if( flag == -2 ) return -2;

	//compute w_mean
	yo -= y_mean;
	mul( yo, C, temp );
	mul( temp, Pw, w_mean );

	//relaxation
	if( inflation_method == 'r' )
	{
		mul( 1. - inflation, W );
		for( int i = 0; i < k; ++i ) W[i][i] += inflation;
	}
	
	//add w_mean to W
	for( int i = 0; i < k; ++i ) W[i] += w_mean;

	//compute analysis ensemble
	mul( W, X, Xa );
	final_Xa:
	for( int i = 0; i < k; ++i ) Xa[i] += x_mean;

	//insert regional analysis ensemble into global ensemble
	for( int i = 0; i < k; ++i )
		for( int j = 0; j < m; ++j )
			Xg[i][ j+region_begin ] = Xa[i][j];

	return 0;
}

inline void ETKF_regional::get_region( const vector &global, vector &regional )
{
	for( int i = 0; i < regional.size(); ++i )
		regional[ i ] = global[ region_begin + i ];
}

// Ensemble transform Kalman filter for joint states with artificial observations of difference between global and regional

ETKF_joint::ETKF_joint( double epsilon_g, double epsilon_r, char inflation_method, int region_begin, int region_end, double Lambda, double Kappa ) : inflation_method( inflation_method ), region_begin( region_begin ), region_end( region_end ), Lambda( Lambda ), Kappa( Kappa )
{
	if( inflation_method == 'x' )
	{
		inflation_g = 1. + epsilon_g;
		inflation_r = 1. + epsilon_r;
	}
	else if( inflation_method == 'r' )
	{
		inflation_g = epsilon_g;
		inflation_r = epsilon_r;
	}
	else { std::cerr << "not a valid inflation method\n"; exit( -1 ); }
}

int ETKF_joint::analysis( matrix &Xg, matrix &Xr, const observation &obsv )
{
	int flag, d = Xr.size_col() / Xg.size_col();
	// the first and the last grid points in the sub gegion in terms of global coordinates
	int mg0 = ( region_begin + d - 1 ) / d, mg1 = region_end / d;

	// dimensions
	int k = Xg.size_row(); // number of ensemble members
	int mg = Xg.size_col(); // dimension of the global state
	int mr = region_end - region_begin + 1; // dimension of the regional state
	int m = mg + mr; // dimension of the combined state
	int l = obsv.size(); // dimension of observation
	int l_total = l + mg1 - mg0 + 1;
	// for ensemble computation
	matrix X( k, m ), Xa( k, m ); // background and analysis ensembles for combined states
	matrix Y( k, l_total ); // observation ensemble corresponding to the state ensemble
	vector x_mean( m, 0. ); // mean of state ensemble
	vector y_mean( l_total, 0. ); // mean of observation ensemble
	vector yg( l ), yr( l ); // temporary global and regional observation
	const int *ob_location = obsv.get_location();
	double r = obsv.covariance();
	
	if( inflation_method == 'x' )
	{
		vector xg_mean( Xg.size_col(), 0. ), xr_mean( Xr.size_col(), 0. );

		for( int i = 0; i < k; ++i )
		{
			xg_mean += Xg[i];
			xr_mean += Xr[i];
		}
		xg_mean /= k;
		xr_mean /= k;
		
		for( int i = 0; i < k; ++i )
		{
			Xg[i] = inflation_g * ( Xg[i]-xg_mean ) + xg_mean;
			Xr[i] = inflation_r * ( Xr[i]-xr_mean ) + xr_mean;
		}
	}
	
	// form a combined background ensemble X combining Xg and Xr
	for( int e=0; e<k; ++e )
	{
		for( int j=0; j<mg; ++j ) X[e][j] = Xg[e][j];
		for( int j=region_begin; j<=region_end; ++j ) X[e][j+mg-region_begin] = Xr[e][j];
	}

	// make Y and compute x_mean and y_mean
	for( int e = 0; e < k; ++e )
	{
		obsv.state_to_ob( Xg[e], yg );
		obsv.state_to_ob( Xr[e], yr );

		for( int o = 0; o < l; ++o )
		{
			if( ob_location[o] >= region_begin && ob_location[o] <= region_end )
				Y[e][o] = (1.-Lambda) * yg[o] + Lambda * yr[o];
			else Y[e][o] = yg[o];
		}
		
		for( int o = l; o < l_total; ++o ) Y[e][o] = Xg[e][ mg0+o-l ] - Xr[e][ d*(mg0+o-l) ];
		
		x_mean += X[e];
		y_mean += Y[e];
	}
	x_mean /= k;
	y_mean /= k;

	// subtract x_mean and y_mean from X and Y
	for( int e = 0; e < k; ++e )
	{
		X[e] -= x_mean;
		Y[e] -= y_mean;
	}
	
	vector yo_temp = obsv.get_observation(), yo( l_total );
	matrix C(k, l_total), Pw( k, k );
	matrix M( k, k ); // for eigenvectors
	matrix Minv( k, k ), W( k, k );
	vector w_mean( k ), wx_mean( m ), temp( k ), e( k ); // for eigenvalues
	
	for( int o = 0; o < l; ++o ) yo[o] = yo_temp[o];
	for( int o = l; o < l_total; ++o ) yo[o] = 0.;
	
	// compute C
	for( int i = 0; i < k; ++i )
	{
		for( int j = 0; j < l; ++j ) C[i][j] = Y[i][j] / r;
		for( int j = l; j < l_total; ++j ) C[i][j] = Y[i][j] * Kappa;
	}	
	
	transpose( C );

	// compute Pw and W
	mul( Y, C, Pw );
	for(int i=0; i<k; ++i) Pw[i][i] += (k-1);

	flag = eigen( Pw, e, M );
	if( flag == -1 ) return -1;

	inverse( M, Minv );
	if( flag == -2 ) return -2;

	for( int i=0; i<k; ++i ) for(int j=0; j<k; ++j) M[i][j] *= sqrt((k-1)/e[j]);
	mul( M, Minv, W );
	inverse( Pw );
	if( flag == -2 ) return -2;

	//compute w_mean
	yo -= y_mean;
	mul( yo, C, temp );
	mul( temp, Pw, w_mean );

	//compute w_mean in model space
	mul( w_mean, X, wx_mean );
	
	//compute analysis perturbations
	mul( W, X, Xa );

	//relaxation
	if( inflation_method == 'r' )
	{
		for( int i = 0; i < k; ++i )
		{
			for( int j = 0; j < mg; ++j )
				Xa[i][j] = ( 1. - inflation_g ) * Xa[i][j] + inflation_g * X[i][j];
			for( int j = mg; j < m; ++j )
				Xa[i][j] = ( 1. - inflation_r ) * Xa[i][j] + inflation_r * X[i][j];
		}
	}
	
	//add wx_mean to Xa
	for( int i=0; i<k; ++i ) Xa[i] += wx_mean;

	for( int i=0; i<k; ++i ) Xa[i] += x_mean;

	// split analysis ensemble Xa into Xg and Xr
	for( int i=0; i<k; ++i )
	{
		for( int j=0; j<mg; ++j ) Xg[i][j] = Xa[i][j];
		for( int j=region_begin; j<=region_end; ++j ) Xr[i][j] = Xa[i][j+mg-region_begin];
	}
	
	return 0;
}

// Local ensemble transform Kalman filter

LETKF::LETKF( int span, int span_avg, double epsilon, char inflation_method )
: span( span ), weights( 0 ), inflation_method( inflation_method )
{
	if( inflation_method == 'm' || inflation_method == 'x' ) inflation = 1. + epsilon;
	else if( inflation_method == 'r' ) inflation = epsilon;
	else { std::cerr << inflation_method << " is not a valid inflation method\n"; exit( -1 ); }
	
	if( span_avg > 0 )
	{
		weights.resize( 2*span + 1 );
		weights = 0.;
	
		for( int i = -span_avg; i <= span_avg; ++i )
			weights[ i+span ] = 1. / ( 2*span_avg+1 );
	}
}

LETKF::LETKF( int span, const vector &weight_avg, double epsilon, char inflation_method )
: span( span ), weights( weight_avg ), inflation_method( inflation_method )
{
	if( inflation_method == 'm' || inflation_method == 'x' ) inflation = 1. + epsilon;
	else if( inflation_method == 'r' ) inflation = epsilon;
	else { std::cerr << inflation_method << " is not a valid inflation method\n"; exit( -1 ); }	
}

inline void LETKF::get_local( int middle, const vector &global, vector &local )
{
	for( int i=0; i<local.size(); ++i ) local[i] = global[ mod(i + middle - span, global.size()) ];
}

int LETKF::analysis( matrix &Xg, const observation &obsv )
{
	int flag;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	
	//for global computation
	matrix Xag( k, mg, 0 ); //global analysis ensemble
	matrix Yg( k ); //global observation ensemble corresponding to global state ensemble
	vector xg_mean( mg, 0 ); //mean of global state ensemble
	vector yg_mean( lg, 0 ); //mean of global observation ensemble

	// multiplicative inflation
	if( inflation_method == 'x' )
	{
		vector xg_mean_temp( mg, 0. );
		for( int i = 0; i < k; ++i ) xg_mean_temp += Xg[i];
		xg_mean_temp /= k;
		
		for( int i = 0; i < k; ++i ) Xg[i] = inflation * ( Xg[i]-xg_mean_temp ) + xg_mean_temp;
	}
	
	//make Yg and compute xg_mean and yg_mean
	for( int i = 0; i < k; ++i )
	{
		obsv.state_to_ob( Xg[i], Yg[i] );
		xg_mean += Xg[i];
		yg_mean += Yg[i];
	}
	xg_mean /= k;
	yg_mean /= k;

	//subtract xg_mean and yg_mean from Xg and Yg
	for( int i=0; i < k; ++i )
	{
		Xg[i] -= xg_mean;
		Yg[i] -= yg_mean;
	}
	
	//for local computation
	matrix X( k, m ); //local background ensemble
	matrix Y( k ); //local observation ensemble
	vector x_mean( m );
	vector y_mean;
	vector yo; //local observation
	matrix C;
	matrix Pw( k, k );
	matrix M( k, k ); //for eigenvectors
	matrix Minv( k, k );
	vector e( k ); //for eigenvalues
	matrix W( k, k );
	vector w_mean( k );
	vector temp( k ), r_inv;
	matrix Xa( k, m ); //local analysis ensemble

	//compute local analysis ensemble for each grid point
	for( int n = 0; n < mg; ++n )
	{
		//get local state and observation ensemble
		for(int i=0; i<k; ++i)
		{
			get_local( n, Xg[i], X[i] );
			obsv.global_to_local( n, span, mg, Yg[i], Y[i] );
		}

		//get local x_mean, y_mean, yo
		get_local( n, xg_mean, x_mean );
		obsv.global_to_local( n, span, mg, yg_mean, y_mean );
		obsv.global_to_local( n, span, mg, obsv.get_observation(), yo );

		if( yo.size() == 0 ) { Xa = X; goto final_Xa; }

		//compute C
		if( weights.size() == 0 )
		{
			C.resize( k, Y.size_col() );
			obsv.cov_inv( n, span, mg, r_inv );
			for( int i = 0; i < k; ++i ) for( int j = 0; j < Y.size_col(); ++j ) C[i][j] = Y[i][j] * r_inv[j];
		}
		else div( Y, obsv.covariance(), C );
		
		transpose( C );

		//compute Pw and W
		mul(Y, C, Pw);
		if( inflation_method == 'm' ) for(int i=0; i<k; ++i) Pw[i][i] += (k-1)/inflation;
		else for(int i=0; i<k; ++i) Pw[i][i] += (k-1);
		flag = eigen(Pw, e, M);
		if( flag == -1 ) return -1;
		flag = inverse(M, Minv);
		if( flag == -2 ) return -2;
		for(int i=0; i<k; ++i) for(int j=0; j<k; ++j) M[i][j] *= sqrt((k-1)/e[j]);
		mul(M, Minv, W);
		flag = inverse(Pw);
		if( flag == -2 ) return -2;

		//compute w_mean
		yo -= y_mean;
		mul(yo, C, temp);
		mul(temp, Pw, w_mean);

		//relaxation
		if( inflation_method == 'r' )
		{
			mul( 1. - inflation, W );
			for( int i = 0; i < k; ++i ) W[i][i] += inflation;
		}
	
		//add w_mean to W
		for(int i=0; i<k; ++i) W[i] += w_mean;

		//compute analysis ensemble
		mul(W, X, Xa);
		final_Xa:
		for(int i=0; i<k; ++i) Xa[i] += x_mean;

		//contribute to the global analysis ensemble at some grid points around this point
		if( weights.size() == 0 )
			for( int i = 0; i < k; ++i ) Xag[i][n] = Xa[i][span];
		else
		{
			double weight;
			for( int j = -span, jg, jl; j <= span; ++j )
			{
				jl = j + span;
				jg = mod( n+j, mg );
				weight = weights[ jl ];
				for( int i = 0; i < k; ++i ) Xag[i][jg] += weight * Xa[i][jl];
			}
		}
	}

	Xg = Xag;

	return 0;
}

int LETKF::get_span()
{
	return span;
}

// Local ensemble transform Kalman filter for a sub-region

LETKF_regional::LETKF_regional
( int span, int span_avg, double epsilon, char inflation_method, int region_begin, int region_end )
: span( span ), weights( 0 ), inflation_method( inflation_method ), region_begin( region_begin ), region_end( region_end )
{
	if( inflation_method == 'm' || inflation_method == 'x' ) inflation = 1. + epsilon;
	else if( inflation_method == 'r' ) inflation = epsilon;
	else { std::cerr << "not a valid inflation method\n"; exit( -1 ); }
	
	if( span_avg > 0 )
	{
		weights.resize( 2*span + 1 );
		weights = 0.;
	
		for( int i = -span_avg; i <= span_avg; ++i )
			weights[ i+span ] = 1. / ( 2*span_avg+1 );
	}
}

int LETKF_regional::analysis( matrix &Xg, const observation &obsv )
{
	int flag;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	double weight;

	//for global computation
	matrix Xag( k, mg, 0 ); //global analysis ensemble
	matrix Yg( k ); //global observation ensemble corresponding to global state ensemble
	vector xg_mean( mg, 0 ); //mean of global state ensemble
	vector yg_mean( lg, 0 ); //mean of global observation ensemble

	// multiplicative inflation
	if( inflation_method == 'x' )
	{
		vector xg_mean_temp( mg, 0. );
		for( int i = 0; i < k; ++i ) xg_mean_temp += Xg[i];
		xg_mean_temp /= k;
		
		for( int i = 0; i < k; ++i ) Xg[i] = inflation * ( Xg[i]-xg_mean_temp ) + xg_mean_temp;
	}
	
	//make Yg and compute xg_mean and yg_mean
	for( int i = 0; i < k; ++i )
	{
		obsv.state_to_ob( Xg[i], Yg[i] );
		xg_mean += Xg[i];
		yg_mean += Yg[i];
	}	
	xg_mean /= k;
	yg_mean /= k;

	//subtract xg_mean and yg_mean from Xg and Yg
	for( int i = 0; i < k; ++i )
	{
		Xg[i] -= xg_mean;
		Yg[i] -= yg_mean;
	}
	
	//for local computation
	int local_begin, local_end;
	matrix X( k, m ); //local background ensemble
	matrix Y( k ); //local observation ensemble
	vector x_mean( m );
	vector y_mean;
	vector yo; //local observation
	matrix C;
	matrix Pw( k, k );
	matrix M( k, k ); //for eigenvectors
	matrix Minv( k, k );
	vector e( k ); //for eigenvalues
	matrix W( k, k );
	vector w_mean( k );
	vector temp( k ), r_inv;
	matrix Xa( k, m ); //local analysis ensemble

	//compute local analysis ensemble for each grid point
	for( int n = region_begin; n <= region_end; ++n )
	{
		local_begin = n - span;
		local_end = n + span;
		
		/*if( local_begin < region_begin ) local_begin = region_begin;
		if( local_end > region_end ) local_end = region_end;
		
		resize the vectors and matrix for local calculation
		X.resize( k, local_end - local_begin + 1 );
		x_mean.resize( local_end - local_begin + 1 );
		Xa.resize( k, local_end - local_begin + 1 );*/

		//get local state and observation ensemble
		for( int i = 0; i < k; ++i )
		{
			get_local( local_begin, local_end, Xg[i], X[i] );
			obsv.global_to_region( local_begin, local_end, Yg[i], Y[i] );
		}

		//get local x_mean, y_mean, yo
		get_local( local_begin, local_end, xg_mean, x_mean );
		obsv.global_to_region( local_begin, local_end, yg_mean, y_mean );
		obsv.global_to_region( local_begin, local_end, obsv.get_observation(), yo );

		if( yo.size() == 0 ) { Xa = X; goto final_Xa; }

		//compute C
		if( weights.size() == 0 )
		{
			C.resize( k, Y.size_col() );
			obsv.cov_inv( n, span, mg, r_inv );
			for( int i = 0; i < k; ++i ) for( int j = 0; j < Y.size_col(); ++j ) C[i][j] = Y[i][j] * r_inv[j];
		}
		else div( Y, obsv.covariance(), C );
		
		transpose(C);

		//compute Pw and W
		mul(Y, C, Pw);
		if( inflation_method == 'm' ) for(int i=0; i<k; ++i) Pw[i][i] += (k-1)/inflation;
		else for(int i=0; i<k; ++i) Pw[i][i] += (k-1);
		flag = eigen(Pw, e, M);
		if( flag == -1 ) return -1;
		flag = inverse(M, Minv);
		if( flag == -2 ) return -2;
		for(int i=0; i<k; ++i) for(int j=0; j<k; ++j) M[i][j] *= sqrt((k-1)/e[j]);
		mul(M, Minv, W);
		flag = inverse(Pw);
		if( flag == -2 ) return -2;

		//compute w_mean
		yo -= y_mean;
		mul(yo, C, temp);
		mul(temp, Pw, w_mean);

		//relaxation
		if( inflation_method == 'r' )
		{
			mul( 1. - inflation, W );
			for( int i = 0; i < k; ++i ) W[i][i] += inflation;
		}
	
		//add w_mean to W
		for(int i=0; i<k; ++i) W[i] += w_mean;

		//compute analysis ensemble
		mul(W, X, Xa);
		final_Xa:
		for(int i=0; i<k; ++i) Xa[i] += x_mean;

		//contribute to the global analysis ensemble at some grid points around this point
		if( weights.size() == 0 )
			for( int i = 0; i < k; ++i ) Xag[i][n] = Xa[i][n-local_begin];
		else
		{
			for( int j = -span, jg, jl; j <= span; ++j )
			{
				jl = j + span;
				jg = mod( n+j, mg );
				weight = weights[ jl ];
				for( int i = 0; i < k; ++i ) Xag[i][jg] += weight * Xa[i][jl];
			}
		}
	}

	//normalize the values near the boundaries.
	if( weights.size() > 0 )
	{
		for( int e = 0; e < k; ++e ) for( int i = 0; i < span; ++i )
		{
			weight = double(2*span+1)/(i+span+1);
			Xag[e][ region_begin + i ] *= weight; 
			Xag[e][ region_end   - i ] *= weight;
		}
	}

	Xg = Xag;

	return 0;
}

inline void LETKF_regional::get_local( int local_begin, int local_end, const vector &global, vector &local )
{
	for( int i=0; i<local.size(); ++i ) local[i] = global[ mod( local_begin + i, global.size() ) ];
}

// Local ensemble transform Kalman filter for joint states with artificial observations of difference between global and regional

LETKF_joint::LETKF_joint( int span, int span_avg, double epsilon_g, double epsilon_r, char inflation_method, int region_begin, int region_end, double Lambda, double Kappa ) : span( span ), span_avg( span_avg ), inflation_method( inflation_method ), region_begin( region_begin ), region_end( region_end ), Lambda( Lambda ), Kappa( Kappa )
{
	if( inflation_method == 'x' )
	{
		inflation_g = 1. + epsilon_g;
		inflation_r = 1. + epsilon_r;
	}
	else { std::cerr << "not a valid inflation method\n"; exit( -1 ); }
}

void LETKF_joint::analysis( matrix &Xg, matrix &Xr, const observation &obsv )
{
	bool whole_region = false;
	if( region_begin == 0 && region_end == Xr.size_col()-1 ) whole_region = true;

	int k = Xg.size_row(), mg = Xg.size_col(), mr = Xr.size_col();
	int d = mr / mg;
	int mg0 = ( region_begin + d - 1 ) / d, mg1 = region_end / d;
	int ls = obsv.size();
	int lc = mg1 - mg0 + 1;
	const int *ob_location = obsv.get_location();
	vector xg_mean( mg, 0. ), xr_mean( mr, 0. ), ys_mean( ls, 0. ), yc_mean( lc, 0. ), yg( ls ), yr( ls );
	matrix Ys( k, ls ), Yc( k, lc ), Xag( k, mg, 0. ), Xar( k, mr, 0. );
	
	// Compute xg_mean, xr_mean
	for( int i = 0; i < k; ++i )
	{
		xg_mean += Xg[i];
		xr_mean += Xr[i];
	}
	xg_mean /= k;
	xr_mean /= k;
	
	// Inflate the background ensemble on the state space
	if( inflation_method == 'x' )
	{
		for( int i = 0; i < k; ++i )
		{
			Xg[i] = inflation_g * ( Xg[i]-xg_mean ) + xg_mean;
			Xr[i] = inflation_r * ( Xr[i]-xr_mean ) + xr_mean;
		}
	}

	// Construct Ys, Yc, and compute ys_mean, yc_mean
	for( int e = 0; e < k; ++e )
	{
		obsv.state_to_ob( Xg[e], yg );
		obsv.state_to_ob( Xr[e], yr );
		
		// State observations
		for( int o = 0; o < ls; ++o )
		{
			if( ob_location[o] >= region_begin && ob_location[o] <= region_end )
				Ys[e][o] = (1.-Lambda) * yg[o] + Lambda * yr[o];
			else Ys[e][o] = yg[o];
		}
		
		// Correction observations of the difference
		for( int o = 0; o < lc; ++o ) Yc[e][o] = Xg[e][ mg0+o ] - Xr[e][ d*(mg0+o) ];

		ys_mean += Ys[e];
		yc_mean += Yc[e];
	}
	ys_mean /= k;
	yc_mean /= k;

	// Subtract xg_mean, xr_mean, ys_mean, yc_mean
	for( int e = 0; e < k; ++e )
	{
		Xg[e] -= xg_mean;
		Xr[e] -= xr_mean;
		Ys[e] -= ys_mean;
		Yc[e] -= yc_mean;
	}

	int ng_begin, ng_end, nr_begin, nr_end, nc_begin, nc_end;
	int n0 = region_begin + mr, n1 = region_end + mr;
	if( whole_region ) { n0 = 0; n1 = 3*mr; }
	double r = obsv.covariance();
	double weight_g = 1. / ( 2*span_avg/d + 1 ), weight_r = 1. / ( 2*span_avg + 1 );
	matrix X( k ), Y( k ), Ystate( k ), C, Pw( k, k ), M( k, k ), Minv( k, k ), W(k, k), Xa( k ), Xa1( k, 2*span + 1 );
	vector x_mean, y_mean, y_mean_state, yo, yo_state, e( k ), w_mean( k ), wx_mean, temp( k ), r_inv;

	for( int n = 0 + mr; n < mr + mr; ++n ) // Use positively shifted indices to simplify coding
	{
		// Do analysis on the grid points either on the global or regional models
		if( ( n < n0 || n > n1 ) && n%d != 0 ) continue;

		// Calculate the begin and end indices of the global and regional local patches and correction area
		ng_begin = ( n - span + d - 1 ) / d;
		ng_end   = ( n + span ) / d;
		nr_begin = n - span;
		nr_end   = n + span;
		if( nr_begin < n0 ) nr_begin = n0;
		if( nr_end   > n1 ) nr_end   = n1;
		if( nr_begin > nr_end ) nr_end = nr_begin - 1;
		nc_begin = ( nr_begin + d - 1 ) / d;
		nc_end   = nr_end / d;
		if( nc_begin > nc_end ) nc_end = nc_begin - 1;

		// Total size = global size + regional size
		X.resize ( k, ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		Xa.resize( k, ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		x_mean.resize ( ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		wx_mean.resize( ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		
		// Get the local joint state and the local state observation ensemble
		for( int i = 0; i < k; ++i )
		{
			get_local( Xg[i], ng_begin, ng_end, X[i], 0 );
			get_local( Xr[i], nr_begin, nr_end, X[i], ng_end - ng_begin + 1 );
			obsv.global_to_region( n-span, n+span, Ys[i], Ystate[i] );			
		}
		
		// Combine the state observation ensemble and the correction observation ensemble
		Y.resize( k, Ystate.size_col() + nc_end - nc_begin + 1 );
		for( int i = 0; i < k; ++i ) for( int j = 0; j < Ystate.size_col(); ++j )
			Y[i][j] = Ystate[i][j];
		for( int i = 0; i < k; ++i )
			get_local( Yc[i], nc_begin - mg - mg0, nc_end - mg - mg0, Y[i], Ystate.size_col() );
			
		// Get local x_mean
		get_local( xg_mean, ng_begin, ng_end, x_mean, 0 );
		get_local( xr_mean, nr_begin, nr_end, x_mean, ng_end - ng_begin + 1 );
	
		// Get local y_mean
		obsv.global_to_region( n-span, n+span, ys_mean, y_mean_state );
		y_mean.resize( y_mean_state.size() + nc_end - nc_begin + 1 );
		for( int i = 0; i < y_mean_state.size(); ++i ) y_mean[i] = y_mean_state[i];
		get_local( yc_mean, nc_begin - mg - mg0, nc_end - mg - mg0, y_mean, y_mean_state.size() );

		// Get local yo assuming the correction observations are 0s
		obsv.global_to_region( n-span, n+span, obsv.get_observation(), yo_state );
		yo.resize( yo_state.size() + nc_end - nc_begin + 1 );
		for( int i = 0; i < yo_state.size(); ++i ) yo[i] = yo_state[i];
		for( int i = yo_state.size(); i < yo.size(); ++i ) yo[i] = 0.;
		if( yo.size() == 0 ) { Xa = X; goto final_Xa; }

		// Compute C
		C.resize ( k, Y.size_col() );
		for( int i = 0; i < k; ++i )
		{
			if( span_avg == -1 )
			{
				obsv.cov_inv( n-mr, span, mr, r_inv );
				for( int j = 0; j < Ystate.size_col(); ++j ) C[i][j] = Y[i][j] * r_inv[j];
			}
			else for( int j = 0; j < Ystate.size_col(); ++j ) C[i][j] = Y[i][j] / r;
			for( int j = Ystate.size_col(); j < Y.size_col(); ++j ) C[i][j] = Y[i][j] * Kappa;
		}
		transpose( C );

		// Compute Pw and W
		mul( Y, C, Pw );
		for( int i = 0; i < k; ++i ) Pw[i][i] += (k-1);
		eigen( Pw, e, M );
		inverse( M, Minv );
		for( int i=0; i<k; ++i ) for( int j=0; j<k; ++j ) M[i][j] *= sqrt((k-1)/e[j]);
		mul( M, Minv, W );
		inverse( Pw );

		// Compute w_mean
		yo -= y_mean;
		mul( yo, C, temp );
		mul( temp, Pw, w_mean );

		// Compute w_mean in model space
		mul( w_mean, X, wx_mean );

		// Compute analysis perturbations
		mul( W, X, Xa );

		// add wx_mean to Xa
		for( int i = 0; i < k; ++i ) Xa[i] += wx_mean;

		final_Xa:
		for( int i = 0; i < k; ++i ) Xa[i] += x_mean;
		
		// Extract regional
		int p = n-span - n0;
		for( int i = 0; i < k; ++i )
		{
			for( int j = ng_end - ng_begin + 1; j < Xa.size_col(); ++j )
			{
				if( p < 0 ) Xa1[i][ j - ng_end + ng_begin - 1 - p ] = Xa[i][j];
				else Xa1[i][ j - ng_end + ng_begin - 1 ] = Xa[i][j];
			}
		}

		// Contribute to the global analysis ensemble at some grid points around this point
		if( n%d == 0 )
		{
			for( int i = n/d - span_avg/d, ig, il; i <= n/d + span_avg/d; ++i )
			{
				il = i - n/d + span/d;
				ig = mod( i, mg );
				for( int e = 0; e < k; ++e ) Xag[e][ig] += weight_g * Xa[e][il];
			}
		}

		// Contribute to the regional analysis ensemble at some grid points around this point
		if( n >= n0 && n <= n1 )
		{
			for( int i = n - span_avg, ig, il; i <= n + span_avg; ++i )
			{
				il = i - n + span;
				ig = mod( i, mr );
				for( int e = 0; e < k; ++e ) Xar[e][ig] += weight_r * Xa1[e][il];
			}
		}
	}

	//normalize the regional values near the boundaries.
	if( !whole_region && span_avg > 0 )
	{
		for( int e = 0; e < k; ++e ) for( int i = 0; i < span_avg; ++i )
		{
			weight_r = double( 2*span_avg+1 ) / ( i+span_avg+1 );
			//weight_r = double( 2*span_avg+1 ) / ( i+span_avg+1 + (span_avg-i)/d );
			Xar[e][ region_begin + i ] *= weight_r;
			Xar[e][ region_end   - i ] *= weight_r;
		}
	}

	Xg = Xag;
	Xr = Xar;
}

/*void LETKF_joint::analysis_cost( matrix &Xg, matrix &Xr, const observation &obsv, double &cost1, double &cost2, double &cost3 )
{
	bool whole_region = false;
	if( region_begin == 0 && region_end == Xr.size_col()-1 ) whole_region = true;

	int k = Xg.size_row(), mg = Xg.size_col(), mr = Xr.size_col(), d = mr / mg;
	int mg0 = ( region_begin + d - 1 ) / d, mg1 = region_end / d;
	int ls = obsv.size(), lc = mg1 - mg0 + 1;
	const int *ob_location = obsv.get_location();
	vector xg_mean( mg, 0. ), xr_mean( mr, 0. ), ys_mean( ls, 0. ), yc_mean( lc, 0. ), yg( ls ), yr( ls );
	matrix Ys( k, ls ), Yc( k, lc ), Xag( k, mg, 0. ), Xar( k, mr, 0. );
	
	// xg_mean, xr_mean
	for( int i = 0; i < k; ++i )
	{
		xg_mean += Xg[i];
		xr_mean += Xr[i];
	}
	xg_mean /= k;
	xr_mean /= k;
	
	// Inflate the background ensemble on the state space
	if( inflation_method == 'x' )
	{
		for( int i = 0; i < k; ++i )
		{
			Xg[i] = inflation_g * ( Xg[i]-xg_mean ) + xg_mean;
			Xr[i] = inflation_r * ( Xr[i]-xr_mean ) + xr_mean;
		}
	}

	// Ys, Yc, ys_mean, yc_mean
	for( int e = 0; e < k; ++e )
	{
		obsv.state_to_ob( Xg[e], yg );
		obsv.state_to_ob( Xr[e], yr );
		
		// State observations
		for( int o = 0; o < ls; ++o )
		{
			if( ob_location[o] >= region_begin && ob_location[o] <= region_end )
				Ys[e][o] = (1.-Lambda) * yg[o] + Lambda * yr[o];
			else Ys[e][o] = yg[o];
		}
		
		// Correction observations of the difference
		for( int o = 0; o < lc; ++o ) Yc[e][o] = Xg[e][ mg0+o ] - Xr[e][ d*(mg0+o) ];

		ys_mean += Ys[e];
		yc_mean += Yc[e];
	}
	ys_mean /= k;
	yc_mean /= k;

	// Subtract xg_mean, xr_mean, ys_mean, yc_mean
	for( int e = 0; e < k; ++e )
	{
		Xg[e] -= xg_mean;
		Xr[e] -= xr_mean;
		Ys[e] -= ys_mean;
		Yc[e] -= yc_mean;
	}

	int ng_begin, ng_end, nr_begin, nr_end, nc_begin, nc_end;
	int n0 = region_begin + mr, n1 = region_end + mr;
	if( whole_region ) { n0 = 0; n1 = 3*mr; }
	double r = obsv.covariance();
	double weight_g = 1. / ( 2*span_avg/d + 1 ), weight_r = 1. / ( 2*span_avg + 1 );
	matrix X( k ), Y( k ), Ystate( k ), C, Pw( k, k ), M( k, k ), Minv( k, k ), W(k, k), Xa( k ), Xa1( k, 2*span + 1 );
	vector x_mean, y_mean, y_mean_state, yo, yo_state, e( k ), w_mean( k ), wx_mean, temp( k ), r_inv;

	// For costs
	int size_g = span/d*2 + 1;
	int size_r = span*2   + 1;
	int size = size_g + size_r;
	int count = 0;
	vector xag_mean( size_g ), xar_mean( size_r ), xa_mean( size );
	vector xag_mean_ext( mg ), xar_mean_ext( mr ), ys_ext( ls ), ys, e_joint( size ), we_mean( size );
	matrix P( size, size), M_joint( size, size ), M_joint_inv( size, size );
	cost1 = cost2 = cost3 = 0.;
	
	for( int n = 0 + mr; n < mr + mr; ++n ) // Use positively shifted indices to simplify coding
	{
		// Do analysis on the grid points either on the global or regional models
		if( ( n < n0 || n > n1 ) && n%d != 0 ) continue;

		// begin and end indices of the global and regional local patches and correction area
		ng_begin = ( n - span + d - 1 ) / d;
		ng_end   = ( n + span ) / d;
		nr_begin = n - span;
		nr_end   = n + span;
		if( nr_begin < n0 ) nr_begin = n0;
		if( nr_end   > n1 ) nr_end   = n1;
		if( nr_begin > nr_end ) nr_end = nr_begin - 1;
		nc_begin = ( nr_begin + d - 1 ) / d;
		nc_end   = nr_end / d;
		if( nc_begin > nc_end ) nc_end = nc_begin - 1;

		// Total size = global size + regional size
		X.resize ( k, ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		Xa.resize( k, ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		x_mean.resize ( ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		wx_mean.resize( ng_end - ng_begin + 1 + nr_end - nr_begin + 1 );
		
		// X (local joint state ensemble), Ystate (local state observation ensemble)
		for( int i = 0; i < k; ++i )
		{
			get_local( Xg[i], ng_begin, ng_end, X[i], 0 );
			get_local( Xr[i], nr_begin, nr_end, X[i], ng_end - ng_begin + 1 );
			obsv.global_to_region( n-span, n+span, Ys[i], Ystate[i] );			
		}
		
		// Y: combine Ystate (state observation ensemble) and the correction observation ensemble
		Y.resize( k, Ystate.size_col() + nc_end - nc_begin + 1 );
		for( int i = 0; i < k; ++i ) for( int j = 0; j < Ystate.size_col(); ++j )
			Y[i][j] = Ystate[i][j];
		for( int i = 0; i < k; ++i )
			get_local( Yc[i], nc_begin - mg - mg0, nc_end - mg - mg0, Y[i], Ystate.size_col() );
			
		// local x_mean
		get_local( xg_mean, ng_begin, ng_end, x_mean, 0 );
		get_local( xr_mean, nr_begin, nr_end, x_mean, ng_end - ng_begin + 1 );
	
		// local y_mean
		obsv.global_to_region( n-span, n+span, ys_mean, y_mean_state );
		y_mean.resize( y_mean_state.size() + nc_end - nc_begin + 1 );
		for( int i = 0; i < y_mean_state.size(); ++i ) y_mean[i] = y_mean_state[i];
		get_local( yc_mean, nc_begin - mg - mg0, nc_end - mg - mg0, y_mean, y_mean_state.size() );

		// local yo assuming the correction observations are 0s
		obsv.global_to_region( n-span, n+span, obsv.get_observation(), yo_state );
		yo.resize( yo_state.size() + nc_end - nc_begin + 1 );
		for( int i = 0; i < yo_state.size(); ++i ) yo[i] = yo_state[i];
		for( int i = yo_state.size(); i < yo.size(); ++i ) yo[i] = 0.;
		if( yo.size() == 0 ) { Xa = X; goto final_Xa; }

		// C
		C.resize ( k, Y.size_col() );
		for( int i = 0; i < k; ++i )
		{
			if( span_avg == -1 )
			{
				obsv.cov_inv( n-mr, span, mr, r_inv );
				for( int j = 0; j < Ystate.size_col(); ++j ) C[i][j] = Y[i][j] * r_inv[j];
			}
			else for( int j = 0; j < Ystate.size_col(); ++j ) C[i][j] = Y[i][j] / r;
			for( int j = Ystate.size_col(); j < Y.size_col(); ++j ) C[i][j] = Y[i][j] * Kappa;
		}
		transpose( C );

		// Pw, W
		mul( Y, C, Pw );
		for( int i = 0; i < k; ++i ) Pw[i][i] += (k-1);
		eigen( Pw, e, M );
		inverse( M, Minv );
		for( int i=0; i<k; ++i ) for( int j=0; j<k; ++j ) M[i][j] *= sqrt((k-1)/e[j]);
		mul( M, Minv, W );
		inverse( Pw );

		// w_mean
		yo -= y_mean;
		mul( yo, C, temp );
		mul( temp, Pw, w_mean );

		// w_mean in model space
		mul( w_mean, X, wx_mean );

		// Analysis perturbations
		mul( W, X, Xa );

		// Add wx_mean to Xa
		for( int i = 0; i < k; ++i ) Xa[i] += wx_mean;

		final_Xa:
		for( int i = 0; i < k; ++i ) Xa[i] += x_mean;
		
		// Extract regional into Xa1 of size 2*span+1. Grid point n corresponds to the middle of Xa1.
		int p = n-span - n0;
		for( int i = 0; i < k; ++i )
		{
			for( int j = ng_end - ng_begin + 1; j < Xa.size_col(); ++j )
			{
				if( p < 0 ) Xa1[i][ j - ng_end + ng_begin - 1 - p ] = Xa[i][j];
				else Xa1[i][ j - ng_end + ng_begin - 1 ] = Xa[i][j];
			}
		}

		// Contribute to the global analysis ensemble at some grid points around this point
		if( n%d == 0 )
		{
			for( int i = n/d - span_avg/d, ig, il; i <= n/d + span_avg/d; ++i )
			{
				il = i - n/d + span/d;
				ig = mod( i, mg );
				for( int e = 0; e < k; ++e ) Xag[e][ig] += weight_g * Xa[e][il];
			}
		}

		// Contribute to the regional analysis ensemble at some grid points around this point
		if( n >= n0 && n <= n1 )
		{
			for( int i = n - span_avg, ig, il; i <= n + span_avg; ++i )
			{
				il = i - n + span;
				ig = mod( i, mr );
				for( int e = 0; e < k; ++e ) Xar[e][ig] += weight_r * Xa1[e][il];
			}
		}
		
		// Compute the 3 terms of the cost function on the global grid points
		// and when the whole regioanl local patch is inside the subregion.
		if( n%d == 0 && nr_begin == n-span && nr_end == n+span )
		{
			// xa_mean (analysis mean), xag_mean (global analysis mean), xar_mean (regional analysis mean)
			xa_mean = wx_mean + x_mean;
			for( int i = 0; i < size_g; ++i ) xag_mean[i] = xa_mean[i];
			for( int i = 0; i < size_r; ++i ) xar_mean[i] = xa_mean[size_g+i];

			// P
			P = 0.;
			for( int i = 0; i < size; ++i ) for( int j = 0; j < size; ++j )
			{
				for( int e = 0; e < k; ++e ) P[i][j] += X[e][i] * X[e][j];
			}
			P /= k-1;

			// cost1: (xa_mean-x_mean)^T Pinv (xa_mean-x_mean)
			eigen( P, e_joint, M_joint );
			bubble_sort( e_joint, M_joint );
			inverse( M_joint, M_joint_inv );
			mul( M_joint_inv, wx_mean, we_mean );
			for( int i = size-1; i > size-1 - (k-1); --i )
			{
				cost1 += we_mean[i] * we_mean[i] / e_joint[i];
			}
			
			// ys = H(xa_mean)
			for( int i = 0; i < size_g; ++i ) xag_mean_ext[ mod(ng_begin+i,mg) ] = xag_mean[i];
			for( int i = 0; i < size_r; ++i ) xar_mean_ext[ mod(nr_begin+i,mr) ] = xar_mean[i];
			obsv.state_to_ob( xag_mean_ext, yg );
			obsv.state_to_ob( xar_mean_ext, yr );
			for( int o = 0; o < ls; ++o ) ys_ext[o] = (1.-Lambda) * yg[o] + Lambda * yr[o];
			obsv.global_to_region( n-span, n+span, ys_ext, ys );

			// cost2: (yo_state-ys)^T (yo_state-ys)
			for( int i = 0; i < ys.size(); ++i )
			{
				cost2 += ( yo_state[i] - ys[i] ) * ( yo_state[i] - ys[i] );
			}

			// cost3: [ Gg(xag_mean) - Gr(xar_mean) ]^T [ Gg(xag_mean) - Gr(xar_mean) ]
			for( int i = -span/d; i <= span/d; ++i )
			{
				cost3 += ( xag_mean[ i + span/d ] - xar_mean[ i*d + span ] )
					   * ( xag_mean[ i + span/d ] - xar_mean[ i*d + span ] );
			}

			++count;
		}
	}

	// Normalize the regional values near the boundaries.
	if( !whole_region && span_avg > 0 )
	{
		for( int e = 0; e < k; ++e ) for( int i = 0; i < span_avg; ++i )
		{
			weight_r = double( 2*span_avg+1 ) / ( i+span_avg+1 );
			//weight_r = double( 2*span_avg+1 ) / ( i+span_avg+1 + (span_avg-i)/d );
			Xar[e][ region_begin + i ] *= weight_r;
			Xar[e][ region_end   - i ] *= weight_r;
		}
	}

	Xg = Xag;
	Xr = Xar;
	
	// Normalize the costs
	cost2 /= r;
	cost3 *= Kappa;
	cost1 /= count;
	cost2 /= count;
	cost3 /= count;
}

void LETKF_joint::bubble_sort( vector &v, matrix &m )
// sort v and corresponding m in an increasing order
{
	bool flag;
	double temp;
	
	do {
		flag = false;
		for(int j=0; j<v.size()-1; ++j)
		{
			if( v[j] > v[j+1] ) 
			{
				temp = v[j]; v[j] = v[j+1]; v[j+1] = temp;
				for(int i=0; i<m.size_row(); ++i)
				{
					temp = m[i][j]; m[i][j] = m[i][j+1]; m[i][j+1] = temp;
				}
				flag = true;
			}
		}
	} while( flag );
}*/

inline void LETKF_joint::get_local( const vector &global, int begin, int end, vector &local, int start_local )
{
	for( int i = begin; i <= end; ++i )	local[ start_local++ ] = global[ mod( i, global.size() ) ];
}

// Local ensemble Kalman filter

LEKF::LEKF( int span, int span_avg, double epsilon ) : span( span ), span_avg( span_avg ), epsilon( epsilon ) {}

inline void LEKF::get_local( int middle, const vector &global, vector &local )
{
	for(int i=0; i<local.size(); ++i) local[i] = global[ mod(i + middle - span, global.size()) ];
}

void LEKF::bubble_sort( vector &v, matrix &m )
// sort v and corresponding m in an increasing order
{
	bool flag;
	double temp;
	
	do {
		flag = false;
		for(int j=0; j<v.size()-1; ++j)
		{
			if( v[j] > v[j+1] ) 
			{
				temp = v[j]; v[j] = v[j+1]; v[j+1] = temp;
				for(int i=0; i<m.size_row(); ++i)
				{
					temp = m[i][j]; m[i][j] = m[i][j+1]; m[i][j+1] = temp;
				}
				flag = true;
			}
		}
	} while( flag );
}

int LEKF::analysis( matrix &Xg_t, const observation &obsv, int k )
// X = sqrt(k')*Xpaper, P = k'*Ppaper, R =k'*Rpaper. The equations remain the same.
{
	int nens = Xg_t.size_row(); // number of ensemble members
	int m = span * 2 +1; // dimension of local state
	if( k==0) k = nens - 1; // estimate of dimenstion of local state
	if( k>m ) k = m;

	int flag;
	double R = (nens-1) * obsv.covariance(), inflation;
	int mg = Xg_t.size_col(); // dimension of global state
	matrix Xag_t(nens, mg, 0); // global analysis ensemble
	
	// xgb_mean : mean of global background ensemble
	vector xgb_mean( mg, 0 );
	for(int i=0; i<nens; ++i) xgb_mean += Xg_t[i];
	xgb_mean /= nens;
	
	// yg_mean : global observation corresponding to xgb_mean
	vector yg_mean;
	obsv.state_to_ob( xgb_mean, yg_mean );

	// Xg_t = Xg_t - xgb_mean
	for(int i=0; i<nens; ++i) Xg_t[i] -= xgb_mean;
	
	// Hg : global Jacobian matrix of observation operator
	matrix Hg(obsv.size(), mg);
	obsv.Jacobian(xgb_mean, Hg);

	// initial Pb : local covariance matrix for the local patch around the last grid point mg-1
	matrix Pb(m, m, 0);
	for(int ig, jg, i=0; i<m; ++i) {
		for(int j=0; j<m; ++j) {
			ig = mod( i + mg-1 - span, mg ); // from local coordinates to the circular global coordinates
			jg = mod( j + mg-1 - span, mg );
			for(int v=0; v<nens; ++v)  Pb[i][j] += Xg_t[v][ig] * Xg_t[v][jg];
		}
	}
	
	// for local computation
	int s; // number of local observations
	vector e_Pb(m); // for the calculation of eigen values of Pb
	matrix U_Pb(m, m); // for the calculation of eigen vectors of Pb
	vector e(nens); // for the calculation of Y
	matrix U(nens, nens), U_i(nens, nens); // for the calculation of Y
	vector lambda(k);
	matrix Q(m, k), Q_t(k, m), Pa_hat(k, k);
	vector xb_mean(m); // local background mean
	vector xa_mean(m); // local analysis mean
	matrix Xb_t(nens, m), Xa_t(nens, m), Xb_hat(k, nens), Xb_hat_t(nens, k), Xa_hat_t(nens, k);
	matrix Y(nens, nens), temp(k, k), temp1(nens, nens);

	// compute local analysis ensemble for each grid point
	for(int n=0; n<mg; ++n)
	{
		vector y; // local observation
		vector y_mean; // local observation corresponding to xb_mean
		matrix H; // local Jacobian matrix
		matrix H_hat; // local Jacobian in eigen basis
		matrix H_hat_t; // transpose

		// Xb_t, xb_mean
		for(int i=0; i<nens; ++i) get_local(n, Xg_t[i], Xb_t[i]);
		get_local(n, xgb_mean, xb_mean);
	
		// y, y_mean
		obsv.global_to_local(n, span, mg, obsv.get_observation(), y);
		obsv.global_to_local(n, span, mg, yg_mean, y_mean);
		if ( (s = y.size()) == 0 ) { Xa_t = Xb_t; goto contribute; }

		// Pb
		// compute Pb from the one at the previous grid point	
		for(int i=0; i<m-1; ++i)
			for(int j=0; j<m-1; ++j)
				Pb[i][j] = Pb[i+1][j+1];
		// compute the remaining column of Pb
		for(int j=m-1, jg = mod( j + n - span, mg ), i=0, ig; i<m-1; ++i) {
			ig = mod( i + n - span, mg );
			Pb[i][j] = 0;
			for(int v=0; v<nens; ++v) Pb[i][j] += Xg_t[v][ig] * Xg_t[v][jg];
		}
		// compute the remaining row of Pb 
		for(int i=m-1, ig = mod( i + n - span, mg ), j=0, jg; j<m; ++j) {
			jg = mod( j + n - span, mg );
			Pb[i][j] = 0;
			for(int v=0; v<nens; ++v) Pb[i][j] += Xg_t[v][ig] * Xg_t[v][jg];
		}

		// Q, Q_t, lambda in a decreasing order
		flag = eigen(Pb, e_Pb, U_Pb);
		if( flag == -1 ) return -1;
		bubble_sort(e_Pb, U_Pb); // sort in an increasing order because it is a typical order of eigen values 
		for(int j=0; j<k; ++j) { // take the largest k eigen values and the corresponding eigen vectors
			lambda[j] = e_Pb[m-1-j];
			for(int i=0; i<m; ++i) Q[i][j] = U_Pb[i][m-1-j];
		}
		transpose( Q, Q_t );

		// Xb_hat, Xb_hat_t : X = Q X_hat, X_hat = Q_t X
		mul( Xb_t, Q, Xb_hat_t );
		transpose( Xb_hat_t, Xb_hat );

		// H_hat and H_hat_t : H_hat = H * Q
		obsv.global_to_local_Jacobian(n, span, mg, Hg, H);
		mul(H, Q, H_hat);
		transpose(H_hat, H_hat_t);
		
		// Pa_hat
		for(int i=0; i<k; ++i) {
			for(int j=0; j<k; ++j) {
				Pa_hat[i][j] = 0;
				for(int v=0; v<s; ++v) Pa_hat[i][j] += H_hat[v][i] * H_hat[v][j];
				Pa_hat[i][j] /= R;
			}
		}
		for(int i=0; i<k; ++i) Pa_hat[i][i] += 1./lambda[i];
		inverse(Pa_hat);
		// enhanced variance inflation
		inflation = epsilon * trace(Pa_hat) / k;
		for(int i=0; i<k; ++i) Pa_hat[i][i] += inflation;
		
		// xa_mean
		xa_mean = Q * Pa_hat * H_hat_t * ( y - y_mean ) / R + xb_mean;

		// Y
		// let diag = diag( 1/lambda1, ..., 1/lambdak )
		// Y^2 = I + Xb_hat_t * diag * (Pa_hat * diag - I ) * Xb_hat 
		for(int i=0; i<k; ++i) for(int j=0; j<k; ++j) temp[i][j] = Pa_hat[i][j] / lambda[j];
		for(int i=0; i<k; ++i) temp[i][i] -= 1.;
		for(int i=0; i<k; ++i) for(int j=0; j<k; ++j) temp[i][j] /= lambda[i];
		Y = Xb_hat_t * temp * Xb_hat;
		for(int i=0; i<nens; ++i) Y[i][i] += 1.; // Y^2
		flag = eigen(Y, e, U);
		if( flag == -1 ) return -1;
		inverse(U, U_i);
		// Y' = diag( sqrt(e_1), ..., sqrt(e_nens) ), Y = U Y' U_i
		for(int i=0; i<nens; ++i) for(int j=0; j<nens; ++j) temp1[i][j] = sqrt(e[i]) * U_i[i][j]; // temp1 = Y' U_i
		mul( U, temp1, Y );

		// Xa_hat_t : = Y * Xb_hat_t
		mul( Y, Xb_hat_t, Xa_hat_t );

		// Xa_t : Xa_t = Xa_hat_t * Q_t;
		mul( Xa_hat_t, Q_t, Xa_t );
		//mul( Y, Xb_t, Xa_t ); // test

		// Add xa_mean to each of the ensemble members of Xa_t
		for(int i=0; i<nens; ++i) Xa_t[i] += xa_mean;
		
		// contribute to the global analysis ensemble at some grid points around this point
		contribute:
		for(int jg, jl, j= -span_avg; j<= span_avg; ++j) {
			jg = mod(n+j, mg);
			jl = span+j;
			for(int i=0; i<nens; ++i) Xag_t[i][jg] += Xa_t[i][jl];
		}
	}

	div(Xag_t, 2*span_avg+1, Xg_t);
	return 0;
}

// Direct interpolation

void interpolation::analysis( vector &x, const observation &obsv )
{
	int n_ob = obsv.size();
	int N = obsv.size_true_state();
	vector ob = obsv.get_observation();
	vector a(N);

	int n0, n1, n2, n3;
	double z0, z1, z2, z3;
	int Lag0, Lag1, Lag2, Lag3;

	for(int i=0; i<n_ob; ++i)
	{
		n0 = obsv.get_location()[i];
		n1 = obsv.get_location()[(i+1)%n_ob];
		n2 = obsv.get_location()[(i+2)%n_ob];
		n3 = obsv.get_location()[(i+3)%n_ob];

		z0 = ob[i];
		z1 = ob[(i+1)%n_ob];
		z2 = ob[(i+2)%n_ob];
		z3 = ob[(i+3)%n_ob];

		n1 = n0 + mod(n1-n0, N);
		n2 = n0 + mod(n2-n0, N);
		n3 = n0 + mod(n3-n0, N);
	
		a[n0] = z0;
		Lag0 = (n0-n1)*(n0-n2)*(n0-n3);
		Lag1 = (n1-n0)*(n1-n2)*(n1-n3);
		Lag2 = (n2-n0)*(n2-n1)*(n2-n3);
		Lag3 = (n3-n0)*(n3-n1)*(n3-n2);

		for(int n=n1+1; n<n2; ++n) {
			a[ mod(n, N) ] = z0*(n-n1)*(n-n2)*(n-n3)/Lag0
				+ z1*(n-n0)*(n-n2)*(n-n3)/Lag1
				+ z2*(n-n0)*(n-n1)*(n-n3)/Lag2
				+ z3*(n-n0)*(n-n1)*(n-n2)/Lag3;
		}
	}

	for(int d = a.size()/x.size(), n=0; n<x.size(); ++n) x[n] = a[n*d];
}
