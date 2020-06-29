#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include "tensor.h"
#include "model.h"
#include "observation.h"
#include "analysis.h"
#define FORECAST
#define RMS_SPATIAL
#define ANALYSIS_LETKF

const char file_out[] = "joint-state-tm3_960_b10_c0.6-gm2_240-Fr14-nRKts48-o15-n1.0-LETKFx-span40-nens40-eg.024-er.02-ba10-K0.01-L0.7-idum-2.tsv";
const char inflation_method = 'x';
const int Nt = 960, Ng = 240, Nr = 960, No = 15, span = 40, nens = 40, nRKts = 48, cycle[] = { 1000, 100000 }, n_forecast = 5;
const int region_begin = /*0*/240, region_end = /*Nr-1*/720, boundary_radius = 10; // in grid points of true model
const double noise = 1., epsilon_g = 0.024, epsilon_r = 0.02, Kappa = 0.01, Lambda = 0.7;
const double F_t = 15., F_r = 14., b_t = 10., b_r = 10., c_t = 0.6, c_r = 0.6, dt = 0.05/nRKts;

void set_ens( const char fine_in[], matrix &ens_g );
void tie( const vector &state_g, vector &state_r, int region_begin, int region_end );
void tie( const LM &model_g, LM &model_r, int region_begin, int region_end );
void rms_error( const matrix &ens_g, const matrix &ens_r, const LM &model_true, double Lambda, double &rms_g, double &rms_r, 	double &rms_sum );
void rms_error( const matrix &ens_g, const matrix &ens_r, const LM &model_true, double Lambda, vector &rms_g, vector &rms_r, 	vector &rms_sum );
void rms_error( const LM &model_g, const LM &model_r, const LM &model_true, vector &rms_g, vector &rms_r );
void process_boundaries( const matrix &ens_g, matrix &ens_r, int region_begin, int region_end, int boundary_radius );
void forecast( const LM &model_true, const matrix &ens_g, const matrix &ens_r, LM &model_t, LM &model_g, LM &model_r, matrix 	&rms_g, matrix &rms_r );
void global_to_region( const vector &state_g, vector &state_r );
void evolve( double time, LM &model_true, LM &model_g, LM &model_r, matrix &ens_g, matrix &ens_r );
void cal_mean( const matrix &ens, vector &mean );
void region_to_global( const matrix &ens_r, matrix &ens_g, int region_begin, int region_end );
void print( const double &rms_g, const double &rms_r, const double &rms_sum, int stage, int c );
void print( const vector &rms_g, const vector &rms_r, const vector &rms_sum, int stage, int c );
void print( const matrix &rms_g, const matrix &rms_r, int stage, int c );
int check_blowup( const matrix &ens_g, const matrix &ens_r, int c );
int check_blowup( const LM &model_g, const LM &model_r, int c );

int main()
{
	int flag, blowup_count = 0;
	std::ofstream fout;
	
	/*fout.open( file_out);
	fout << "# global\tregional\tcombined\n";
	fout.close();*/
	/*fout.open( file_out );
	fout << "# sum of the sqaures of the errors at every 1000 cycles descarding the first 1000 cycles\n"
		 << "# cycle_id\tsum-global\t\tsum-regional\t\tsum-combined\n";
	fout.close();*/
	
	#ifdef FORECAST
	matrix rms_g( n_forecast, Ng, 0. ), rms_r( n_forecast, Nr, 0. ), rms_sum( n_forecast, Nr, 0. );
	#else
		#ifdef RMS_SPATIAL
		vector rms_g( Ng, 0. ), rms_r( Nr, 0. ), rms_sum( Nr, 0. );
		#else
		double rms_g = 0., rms_r = 0., rms_sum = 0.;
		#endif
	#endif
	std::vector<int> loc( No );
	choose_regular_numbers( No, Nt, loc );
	LM model_true( 3, Nt, F_t, dt, b_t, c_t ), model_t( 3, Nt, F_t, dt, b_t, c_t );
	LM model_r( 3, Nr, F_r, dt, b_r, c_r ), model_g( 2, Ng, F_t, dt );
	observation obsv( noise, loc, Nt );
	#if defined ANALYSIS_ETKF
	ETKF_joint analysis( epsilon_g, epsilon_r, inflation_method, region_begin, region_end, Lambda, Kappa );
	#elif defined ANALYSIS_LETKF
	LETKF_joint analysis( span, 0, epsilon_g, epsilon_r, inflation_method, region_begin, region_end, Lambda, Kappa );
	#endif
	matrix ens_g( nens, Ng ), ens_r( nens, Nr );

	std::ifstream fin( "ensemble-preparation-m3_960-nens960" );
	for( int i = 0; i < 2; ++i ) for( int n = 0; n < Nt; ++n ) fin >> model_true[n];
	fin.close();

	set_ens( "ensemble-preparation-m2_240-nens480", ens_g );
	for( int e = 0; e < nens; ++e ) global_to_region( ens_g[e], ens_r[e] );
	
	for( int stage = 0; stage < 2; ++stage, rms_g = 0., rms_r = 0., rms_sum = 0. )
	for( int c = 1; c <= cycle[ stage ]; ++c )
	{
		obsv.make_observation( model_true.get_state() );

		analysis.analysis( ens_g, ens_r, obsv );
		
		if( !( region_begin == 0 && region_end == Nr-1 ) )
			process_boundaries( ens_g, ens_r, region_begin, region_end, boundary_radius );
		
		#ifdef FORECAST
		forecast( model_true, ens_g, ens_r, model_t, model_g, model_r, rms_g, rms_r );
		blowup_count += check_blowup( model_g, model_r, c );
		#else
		rms_error( ens_g, ens_r, model_true, Lambda, rms_g, rms_r, rms_sum );
		#endif
		
		// Print time series of the rms errors
		/*{
			int d = Nt / Ng;
			int Ns = region_end - region_begin + 1;
			std::ofstream fout;
	
			fout.open( file_out, std::ios_base::app );

			fout << sqrt( rms_g / Ng ) << '\t' << sqrt( rms_r / Ns ) << '\t' << sqrt( rms_sum / Ns ) << '\n';
		 
			fout.close();
			
			rms_g = rms_r = rms_sum = 0.;
		}*/
		
		for( int i = 0; i < nRKts; ++i )
		{
			if( !( region_begin == 0 && region_end == Nr-1 ) )
				for( int e = 0; e < nens; ++e ) tie( ens_g[e], ens_r[e], region_begin, region_end );
			evolve( dt, model_true, model_g, model_r, ens_g, ens_r );
		}

		flag = check_blowup( ens_g, ens_r, c );
		if( flag ) return 0;
		
		if( c%1000 == 0 && stage == 1 )
		{
			print( rms_g, rms_r, stage, c );
			fout.open( file_out, std::ios_base::app );
			fout << "# The forecast blew up " << blowup_count << " times\n";
			fout.close();
			//rms_g = rms_r = rms_sum = 0.;
		}
	}
}

void set_ens( const char fine_in[], matrix &ens_g )
{
	std::ifstream fin( fine_in );
	
	for( int e = 0; e < nens; ++e ) for( int n = 0; n < Ng; ++n )
		fin >> ens_g[e][n];
		
	fin.close();
}

int check_blowup( const matrix &ens_g, const matrix &ens_r, int c )
{
	int flag_g = 0, flag_r = 0;
	std::ofstream fout;
	
	for( int e = 0; e < nens; ++e ) for( int n = 0; n < Ng; n += 10 )
		if( ens_g[e][n] != ens_g[e][n] ) flag_g = 1;
	for( int e = 0; e < nens; ++e ) for( int n = 0; n < Nr; n += 10 )
		if( ens_r[e][n] != ens_r[e][n] ) flag_r = 1;
			
	if( flag_g )
	{
		fout.open( file_out, std::ios_base::app );
		fout << "# The global ensemble BLEW UP at " << c << "th cycle.\n\n";
		fout.close();
	}
	
	if( flag_r )
	{
		fout.open( file_out, std::ios_base::app );
		fout << "# The regional ensemble BLEW UP at " << c << "th cycle.\n\n";
		fout.close();
	}

	if( flag_g || flag_r ) return 1;
	else return 0;
}

int check_blowup( const LM &model_g, const LM &model_r, int c )
{
	int flag_g = 0, flag_r = 0;
	std::ofstream fout;
	
	for( int n = 0; n < Ng; n += 10 ) if( model_g[n] != model_g[n] ) { flag_g = 1; break; }
	for( int n = 0; n < Nr; n += 10 ) if( model_r[n] != model_r[n] ) { flag_r = 1; break; }
		
	if( flag_g || flag_r ) return 1;
	else return 0;
}

void tie( const vector &state_g, vector &state_r, int region_begin, int region_end )
{
	int M = Nr / Ng;
	double y0, y1, slope;

	for( int n = 0; n < Ng; ++n )
	{
		y0 = state_g[ n ];
		y1 = n+1 < Ng ? state_g[ n+1 ] : state_g[ 0 ];
		slope = ( y1 - y0 ) / M;
		
		for( int i, m = 0; m < M; ++m )
		{
			i = n*M + m;
			if( i < region_begin || i > region_end ) state_r[ i ] = y0 + slope * m;
		}
	}
}

void tie( const LM &model_g, LM &model_r, int region_begin, int region_end )
{
	int M = Nr / Ng;
	double y0, y1, slope;

	for( int n = 0; n < Ng; ++n )
	{
		y0 = model_g[ n ];
		y1 = n+1 < Ng ? model_g[ n+1 ] : model_g[ 0 ];
		slope = ( y1 - y0 ) / M;
		
		for( int i, m = 0; m < M; ++m )
		{
			i = n*M + m;
			if( i < region_begin || i > region_end ) model_r[ i ] = y0 + slope * m;
		}
	}
}

void region_to_global( const matrix &ens_r, matrix &ens_g, int region_begin, int region_end )
{
	int q = Nr / Ng;

	for( int e = 0; e < nens; ++e )
		for( int m = ( region_begin + q - 1 ) / q; m <= region_end / q; ++m )
			ens_g[e][m] = ens_r[e][m*q];
}

void rms_error( const matrix &ens_g, const matrix &ens_r, const LM &model_true, double Lambda, double &rms_g, double &rms_r, double &rms_sum )
{
	int d = Nt / Ng;
	vector mean_g( Ng, 0. ), mean_r( Nr, 0. ), mean_g2r( Nr ), mean_sum( Nr );

	for( int e = 0; e < nens; ++e ) { mean_g += ens_g[e]; mean_r += ens_r[e]; }
	mean_g /= nens; mean_r /= nens;
	
	global_to_region( mean_g, mean_g2r );
	mean_sum = Lambda * mean_r + ( 1. - Lambda ) * mean_g2r;
	
	for( int n = 0; n < Ng; ++n )
		rms_g += ( mean_g[n] - model_true[n*d] ) * ( mean_g[n] - model_true[n*d] );
	
	for( int n = region_begin; n <= region_end; ++n )
		rms_r += ( mean_r[n] - model_true[n] ) * ( mean_r[n] - model_true[n] );

	for( int n = region_begin; n <= region_end; ++n )
		rms_sum += ( mean_sum[n] - model_true[n] ) * ( mean_sum[n] - model_true[n] );
}

void rms_error( const matrix &ens_g, const matrix &ens_r, const LM &model_true, double Lambda, vector &rms_g, vector &rms_r, vector &rms_sum )
{
	int d = Nt / Ng;
	vector mean_g( Ng, 0. ), mean_r( Nr, 0. ), mean_g2r( Nr ), mean_sum( Nr );

	for( int e = 0; e < nens; ++e ) { mean_g += ens_g[e]; mean_r += ens_r[e]; }
	mean_g /= nens; mean_r /= nens;
	
	global_to_region( mean_g, mean_g2r );
	mean_sum = Lambda * mean_r + ( 1. - Lambda ) * mean_g2r;
	
	for( int n = 0; n < Ng; ++n )
		rms_g[n] += ( mean_g[n] - model_true[n*d] ) * ( mean_g[n] - model_true[n*d] );
	
	for( int n = region_begin; n <= region_end; ++n )
		rms_r[n] += ( mean_r[n] - model_true[n] ) * ( mean_r[n] - model_true[n] );

	for( int n = region_begin; n <= region_end; ++n )
		rms_sum[n] += ( mean_sum[n] - model_true[n] ) * ( mean_sum[n] - model_true[n] );
}

void rms_error( const LM &model_g, const LM &model_r, const LM &model_true, vector &rms_g, vector &rms_r )
{
	int d = Nt / Ng;
		
	for( int n = 0; n < Ng; ++n )
		rms_g[n] += ( model_g[n] - model_true[n*d] ) * ( model_g[n] - model_true[n*d] );
	
	for( int n = region_begin; n <= region_end; ++n )
		if( model_r[n] > -40 && model_r[n] < 40 )
			rms_r[n] += ( model_r[n] - model_true[n] ) * ( model_r[n] - model_true[n] );
		else rms_r[n] += 10.;
}

void global_to_region( const vector &state_g, vector &state_r )
{
	int M = Nr / Ng;
	double y0, y1, slope;

	for( int n = 0; n < Ng; ++n )
	{
		y0 = state_g[ n ];
		y1 = n+1 < Ng ? state_g[ n+1 ] : state_g[ 0 ];
		slope = ( y1 - y0 ) / M;
		
		for( int i, m = 0; m < M; ++m )
		{
			i = n*M + m;
			state_r[ i ] = y0 + slope * m;
		}
	}
}

void process_boundaries( const matrix &ens_g, matrix &ens_r, int region_begin, int region_end, int boundary_radius )
{
	vector temp( Nr );
	double c;
	
	for( int e = 0; e < nens; ++e )
	{
		global_to_region( ens_g[e], temp );
		
		for( int n = 0; n < boundary_radius; ++n )
		{
			c = double(n) / boundary_radius;
			ens_r[e][ region_begin + n ] = c * ens_r[e][ region_begin + n ] + (1.-c) * temp[ region_begin + n ];
			ens_r[e][ region_end - n ] = c * ens_r[e][ region_end - n ] + (1.-c) * temp[ region_end - n ];
		}
	}
}

void evolve( double time, LM &model_true, LM &model_g, LM &model_r, matrix &ens_g, matrix &ens_r )
{
	model_true.evolve( time );
	
	for( int i = 0; i < nens; ++i )
	{
		model_g.set_state( ens_g[i] );
		model_g.evolve( time );
		model_g.get_state( ens_g[i] );
	}

	for( int i = 0; i < nens; ++i )
	{
		model_r.set_state( ens_r[i] );
		model_r.evolve( time );
		model_r.get_state( ens_r[i] );
	}
}

void forecast( const LM &model_true, const matrix &ens_g, const matrix &ens_r, LM &model_t, LM &model_g, LM &model_r, matrix &rms_g, matrix &rms_r )
{
	vector mean_g( Ng ), mean_r( Nr ), temp_t( Nt ), mean_g2r( Nr );
	
	cal_mean( ens_g, mean_g );
	cal_mean( ens_r, mean_r );
	
	global_to_region( mean_g, mean_g2r );
	mean_r = Lambda * mean_r + ( 1. - Lambda ) * mean_g2r;

	model_true.get_state( temp_t );
	model_t.set_state( temp_t );
	model_g.set_state( mean_g );
	model_r.set_state( mean_r );
	
	rms_error( model_g, model_r, model_t, rms_g[0], rms_r[0] );

	for( int t = 1; t < n_forecast; ++t )
	{
		for( int i = 0; i < nRKts; ++i )
		{
			model_t.evolve( dt );
			if( !( region_begin == 0 && region_end == Nr-1 ) )
				tie( model_g, model_r, region_begin, region_end );
			model_g.evolve( dt );
			model_r.evolve( dt );
		}
		
		rms_error( model_g, model_r, model_t, rms_g[t], rms_r[t] );
	}
}

void cal_mean( const matrix &ens, vector &mean )
{
	int nrow = ens.size_row(), ncol = ens.size_col();

	mean = 0.;
	for( int r = 0; r < nrow; ++r ) mean += ens[r];
	mean /= nrow;
}

void print( const double &rms_g, const double &rms_r, const double &rms_sum, int stage, int c )
{
	int Ns = region_end - region_begin + 1;
	std::ofstream fout;
	
	//fout.open( file_out );
	fout.open( file_out, std::ios_base::app );
	
	fout.precision( 16 );

	if( stage == 0 ) fout << "# pre : ";

	//fout << "number of cycle : " << c << '\n'
	//	 << "Lambda\tKappa\n"
	//	 << Lambda << '\t' << Kappa << '\n'
	//	 << "sum-global\t\tsum-regional\t\tsum-combined\n"
	//	 << rms_g / Ng << '\t' << rms_r / Ns << '\t' << rms_sum / Ns << "\n\n";
	
	fout << c/1000 << '\t' << rms_g / Ng << '\t' << rms_r / Ns << '\t' << rms_sum / Ns << '\n';
	
	fout.close();
}

void print( const vector &rms_g, const vector &rms_r, const vector &rms_sum, int stage, int c )
{
	int d = Nt / Ng;
	std::ofstream fout( file_out );
	
	if( stage == 0 ) fout << "# pre : ";

	fout << "# number of cycle : " << c << '\n'
		 << "# n\trms-global\trms-regional\trms-combined\n";

	for( int n = 0; n < Nt; ++n )
	{
		fout << n << '\t';
		
		if( n%d == 0 )
		{
			if( n >= region_begin && n <= region_end )
				fout << sqrt( rms_g[n/d]/c ) << '\t' << sqrt( rms_r[n]/c ) << '\t' << sqrt( rms_sum[n]/c ) << '\n';
			else
				fout << sqrt( rms_g[n/d]/c ) << '\t' << "NA\tNA\n";
		}
		else
		{
			if( n >= region_begin && n <= region_end )
				fout << "NA\t" << sqrt( rms_r[n]/c ) << '\t' << sqrt( rms_sum[n]/c ) << '\n';
			else fout << "NA\tNA\tNA\n";
		}
	}
		
	fout.close();
}

void print( const matrix &rms_g, const matrix &rms_r, int stage, int c )
{
	int d = Nt / Ng;//, n = c/1000;
	std::ofstream fout;
	//std::ostringstream sout;
	
	//if( n < 10 ) sout << file_out << "-0" << n << ".tsv";
	//else sout << file_out << "-" << n << ".tsv";
	
	//fout.open( sout.str().c_str() );
	fout.open( file_out );
	
	if( stage == 0 ) fout << "# pre : ";

	fout << "# number of cycle : " << c << '\n'
		 << "# n\trms-global\trms-combined\n";
	
	for( int n = 0; n < Nt; ++n )
	{
		fout << n << '\t';
		
		for( int t = 0; t < n_forecast; ++t )
		{
			if( n%d == 0 )
			{
				if( n >= region_begin && n <= region_end )
					fout << sqrt( rms_g[t][n/d]/c ) << '\t' << sqrt( rms_r[t][n]/c ) << '\t';
				else fout << sqrt( rms_g[t][n/d]/c ) << '\t' << "NA\t";
			}
			else
			{
				if( n >= region_begin && n <= region_end )
					fout << "NA\t" << sqrt( rms_r[t][n]/c ) << '\t';
				else fout << "NA\tNA\t";
			}
		}
		
		fout << '\n';
	}

	fout.close();
}
