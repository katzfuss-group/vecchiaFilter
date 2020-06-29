#include <cmath>
#include <iostream>
#include <fstream>
#include "analysis_mpi.h"

const int data_size = 10;

// Local Ensemble Transform Kalman Filter

LETKF_master::LETKF_master( int span, int span_avg, double epsilon, int np, int rank ) 
: span(span), span_avg(span_avg), rho(1+epsilon), weights(0), np(np), rank(rank)
{
	set_weights();
}

LETKF_worker::LETKF_worker( int span, int span_avg, double epsilon, int np, int rank ) 
: span(span), span_avg(span_avg), rho(1+epsilon), weights(0), np(np), rank(rank)
{
	set_weights();
}

void LETKF_master::set_parameters(int span, int span_avg, double epsilon )
{
	LETKF_master::span = span;
	LETKF_master::span_avg = span_avg;
	rho = 1 + epsilon;
	set_weights();
}

void LETKF_worker::set_parameters(int span, int span_avg, double epsilon )
{
	LETKF_worker::span = span;
	LETKF_worker::span_avg = span_avg;
	rho = 1 + epsilon;
	set_weights();
}

void LETKF_master::set_weights()
{
	weights.resize( 2*span_avg + 1 );
	
	/*int weight_total = 0;
	for( int temp, i = -span_avg; i <= span_avg; ++i ) {
		temp = span_avg+1-abs(i);
		weights[i+span_avg] = temp;
		weight_total += temp;
	}
	
	weights /= weight_total;*/
	weights = 1. / weights.size();
	/*weights = 0;
	for( int i=0; i<=span_avg; ++i )
		weights[span_avg+i] = 1./(span_avg+1);*/
}

void LETKF_worker::set_weights()
{
	weights.resize( 2*span_avg + 1 );
	
	/*int weight_total = 0;
	for( int temp, i = -span_avg; i <= span_avg; ++i ) {
		temp = span_avg+1-abs(i);
		weights[i+span_avg] = temp;
		weight_total += temp;
	}
	
	weights /= weight_total;*/
	weights = 1. / weights.size();
	/*weights = 0;
	for( int i=0; i<=span_avg; ++i )
		weights[span_avg+i] = 1./(span_avg+1);*/
}

int LETKF_master::get_span()
{
	return span;
}

int LETKF_master::analysis_ensemble_old(matrix &Xg, const observation &obsv)
{
	int flag = 0;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	MPI::Status status;
	
	//for global computation
	matrix Xag(k, mg, 0); //global analysis ensemble
	matrix Yg(k); //global observation ensemble corresponding to global state ensemble
	vector xg_mean(mg, 0); //mean of global state ensemble
	vector yg_mean(lg, 0); //mean of global observation ensemble
	
	//make Yg and compute xg_mean and yg_mean
	for(int i=0; i<k; ++i)
	{
		obsv.state_to_ob(Xg[i], Yg[i]);
		xg_mean += Xg[i];
		yg_mean += Yg[i];
	}	
	xg_mean /= k;
	yg_mean /= k;

	//subtract xg_mean and yg_mean from Xg and Yg
	for(int i=0; i<k; ++i) {
		Xg[i] -= xg_mean;
		Yg[i] -= yg_mean;
	}

	for(int i=0; i<k; ++i) {
		MPI::COMM_WORLD.Bcast(Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0);
	}
	MPI::COMM_WORLD.Bcast(xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0);

	int n_data = mg / data_size;
	int n_init_send = n_data>np-1 ? np-1 : n_data;
	int n_sent = 0;
	
	// send min(n_data, np-1) data initially
	for(int i=0; i<n_init_send; ++i, ++n_sent)
		MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, i+1, i );

	int buffer_ncols = data_size + 2*span_avg;
	int buffer_size = k * buffer_ncols;
	double *buffer = new double[buffer_size];

	// receive and send more points 
	for(int sender, data_id, rec=0; rec<n_data; ++rec) // receive n_data times
	{
		MPI::COMM_WORLD.Recv(buffer, buffer_size, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
		sender = status.Get_source();
		
		if( n_sent < n_data ) {
			MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_sent++ );
		} else MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_data );

		data_id = status.Get_tag();
		if( data_id > n_data ) { flag = -1; continue; }
		
		//contribute to the global analysis ensemble		
		for(int count=0, jg, j=0; j<buffer_ncols; ++j) {
			jg = cycle( data_id * data_size - span_avg + j, mg );
			for(int i=0; i<k; ++i ) Xag[i][jg] += buffer[count++];
		}
	}
		
	Xg = Xag;
	
	delete [] buffer;
	
	return flag;
}

void LETKF_worker::analysis_ensemble_old(matrix &Xg, const observation &obsv)
{
	int flag;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	int n_data = mg / data_size;
	
	matrix Yg(k, lg);
	vector xg_mean(mg), yg_mean(lg);
	
	//for local computation
	matrix X(k, m); //local background ensemble
	matrix Y(k); //local observation ensemble
	vector x_mean(m);
	vector y_mean;
	vector yo; //local observation
	matrix C;
	matrix Pw(k, k);
	matrix M(k, k); //for eigenvectors
	matrix Minv(k, k);
	vector e(k); //for eigenvalues
	matrix W(k, k);
	vector w_mean(k);
	vector temp(k);
	matrix Xa(k, m); //local analysis ensemble
	
	for(int i=0; i<k; ++i) {
		MPI::COMM_WORLD.Bcast(Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0);
	}
	MPI::COMM_WORLD.Bcast(xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0);
	
	if (rank>n_data ) return;
	
	int buffer_ncols = data_size + 2*span_avg;
	int buffer_size = k * buffer_ncols;
	double *buffer = new double[buffer_size];
	MPI::Status status;

	//compute local analysis ensemble for the partition
	int x, data_id, start_n, end_n;
	while(1)
	{
		MPI::COMM_WORLD.Recv(&x, 0, MPI::INT, 0, MPI::ANY_TAG, status);
		data_id = status.Get_tag();
		if( data_id==n_data ) break;
		
		for(int i=0; i<buffer_size; ++i) buffer[i] = 0;
		int start_n = data_id * data_size;
		int end_n = start_n + data_size;
		for(int buffer_location_start=0, n=start_n; n<end_n; ++n, buffer_location_start += k)
		{
			//get local state and observation ensemble
			for(int i=0; i<k; ++i)
			{
				get_local(n, Xg[i], X[i]);
				obsv.global_to_local(n, span, mg, Yg[i], Y[i]);
			}
			
			//get local x_mean, y_mean, yo
			get_local(n, xg_mean, x_mean);
			obsv.global_to_local(n, span, mg, yg_mean, y_mean);
			obsv.global_to_local(n, span, mg, obsv.get_observation(), yo);

			if (Y.size_col()==0) { Xa = X; goto contribute; } 

			//compute C
			div(Y, obsv.covariance(), C);
			transpose(C);

			//compute Pw and W
			mul(Y, C, Pw);
			for(int i=0; i<k; ++i) Pw[i][i] += (k-1)/rho;
			flag = eigen( Pw, e, M );
			if( flag != 0 ) { data_id = n_data + 1; break; }
			flag = inverse(M, Minv);
			if( flag != 0 ) { data_id = n_data + 1; break; }
			for(int i=0; i<k; ++i) for(int j=0; j<k; ++j) M[i][j] *= sqrt((k-1)/e[j]);
			mul(M, Minv, W);
			flag = inverse(Pw);
			if( flag != 0 ) { data_id = n_data + 1; break; }
		
			//compute w_mean
			yo -= y_mean;
			mul(yo, C, temp);
			mul(temp, Pw, w_mean);

			//add w_mean to W
			for(int i=0; i<k; ++i) W[i] += w_mean;

			//compute analysis ensemble
			mul(W, X, Xa);
			
			contribute:
			for(int i=0; i<k; ++i) Xa[i] += x_mean;
		
			double weight;
			int avg_size = 2*span_avg+1;
			for(int buffer_location=buffer_location_start, jl, j=0; j< avg_size; ++j) {
				weight = weights[j];
				jl = span - span_avg + j;
				for(int i=0; i<k; ++i)  buffer[buffer_location++] += weight * Xa[i][jl];
			}
		}
			
		// send the result
		MPI::COMM_WORLD.Send( buffer, buffer_size, MPI::DOUBLE, 0, data_id );
	}
	
	delete [] buffer;
}

int LETKF_master::analysis_ensemble(matrix &Xg, const observation &obsv)
{
	int flag = 0;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	MPI::Status status;
	
	//for global computation
	matrix Xag(k, mg, 0); //global analysis ensemble
	matrix Yg(k); //global observation ensemble corresponding to global state ensemble
	vector xg_mean(mg, 0); //mean of global state ensemble
	vector yg_mean(lg, 0); //mean of global observation ensemble
	
	//make Yg and compute xg_mean and yg_mean
	for(int i=0; i<k; ++i)
	{
		obsv.state_to_ob(Xg[i], Yg[i]);
		xg_mean += Xg[i];
		yg_mean += Yg[i];
	}	
	xg_mean /= k;
	yg_mean /= k;

	//subtract xg_mean and yg_mean from Xg and Yg
	for(int i=0; i<k; ++i)
	{
		Xg[i] -= xg_mean;
		Yg[i] -= yg_mean;
	}

	for(int i=0; i<k; ++i)
	{
		MPI::COMM_WORLD.Bcast(Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0);
	}
	MPI::COMM_WORLD.Bcast(xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0);

	int n_data = mg;
	int n_init_send = n_data>np-1 ? np-1 : n_data;
	int n_sent = 0;
	
	// send min(n_data, np-1) data initially
	for(int i=0; i<n_init_send; ++i, ++n_sent)
		MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, i+1, i );

	int buffer_ncols = 2 * span_avg + 1;
	int buffer_size = k * buffer_ncols;
	double *buffer = new double[buffer_size];

	// receive and send more points
	for(int sender, data_id, rec=0; rec<n_data; ++rec) // receive n_data times
	{
		MPI::COMM_WORLD.Recv(buffer, buffer_size, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
		sender = status.Get_source();
		
		if( n_sent < n_data ) {
			MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_sent++ );
		} else MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_data );

		data_id = status.Get_tag();
		if( data_id > n_data ) { flag = -1; continue; }
		
		//contribute to the global analysis ensemble
		for(int count=0, jg, j=0; j<buffer_ncols; ++j)
		{
			jg = cycle( data_id - span_avg + j, mg );
			for(int i=0; i<k; ++i ) Xag[i][jg] += buffer[count++];
		}
	}
		
	Xg = Xag;
	
	delete [] buffer;
	
	return flag;
}

void LETKF_worker::analysis_ensemble(matrix &Xg, const observation &obsv )
{
	int flag;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	int n_data = mg;
	
	matrix Yg(k, lg);
	vector xg_mean(mg), yg_mean(lg);
	
	//for local computation
	matrix X(k, m); //local background ensemble
	matrix Y(k); //local observation ensemble
	vector x_mean(m);
	vector y_mean;
	vector yo; //local observation
	matrix C;
	matrix Pw(k, k);
	matrix M(k, k); //for eigenvectors
	matrix Minv(k, k);
	vector e(k); //for eigenvalues
	matrix W(k, k);
	vector w_mean(k);
	vector temp(k);
	matrix Xa(k, m); //local analysis ensemble
	
	for(int i=0; i<k; ++i)
	{
		MPI::COMM_WORLD.Bcast(Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0);
	}
	MPI::COMM_WORLD.Bcast(xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0);
	
	if (rank>n_data ) return;
	
	int buffer_ncols = 2 * span_avg + 1;
	int buffer_size = k * buffer_ncols;
	double *buffer = new double[buffer_size];
	MPI::Status status;

	//compute local analysis ensemble
	int x, data_id, n;
	while(1)
	{
		MPI::COMM_WORLD.Recv(&x, 0, MPI::INT, 0, MPI::ANY_TAG, status);
		data_id = status.Get_tag();
		if( data_id==n_data ) break;
		n = data_id;
		
		//get local state and observation ensemble
		for(int i=0; i<k; ++i)
		{
			get_local(n, Xg[i], X[i]);
			obsv.global_to_local(n, span, mg, Yg[i], Y[i]);
		}
			
		//get local x_mean, y_mean, yo
		get_local(n, xg_mean, x_mean);
		obsv.global_to_local(n, span, mg, yg_mean, y_mean);
		obsv.global_to_local(n, span, mg, obsv.get_observation(), yo);

		if (Y.size_col()==0) { Xa = X; goto contribute; } 

		//compute C
		div(Y, obsv.covariance(), C);
		transpose(C);

		//compute Pw and W
		mul(Y, C, Pw);
		for(int i=0; i<k; ++i) Pw[i][i] += (k-1)/rho;
		flag = eigen( Pw, e, M );
		if( flag != 0 ) { data_id = n_data + 1; break; }
		flag = inverse(M, Minv);
		if( flag != 0 ) { data_id = n_data + 1; break; }
		for(int i=0; i<k; ++i) for(int j=0; j<k; ++j) M[i][j] *= sqrt((k-1)/e[j]);
		mul(M, Minv, W);
		flag = inverse(Pw);
		if( flag != 0 ) { data_id = n_data + 1; break; }
		
		//compute w_mean
		yo -= y_mean;
		mul(yo, C, temp);
		mul(temp, Pw, w_mean);

		//add w_mean to W
		for(int i=0; i<k; ++i) W[i] += w_mean;

		//compute analysis ensemble
		mul(W, X, Xa);
			
		contribute:
		for(int i=0; i<k; ++i) Xa[i] += x_mean;

		double weight;
		int avg_size = 2*span_avg+1, buffer_location = 0;
		for(int jl, j=0; j< avg_size; ++j)
		{
			weight = weights[j];
			jl = span - span_avg + j;
			for(int i=0; i<k; ++i) buffer[buffer_location++] = weight * Xa[i][jl];
		}
			
		// send the result
		MPI::COMM_WORLD.Send( buffer, buffer_size, MPI::DOUBLE, 0, data_id );
	}
	
	delete [] buffer;
}

inline int LETKF_master::cycle(int i, int N)
{
	int cyclic_i = i % N;
	return cyclic_i<0 ? cyclic_i + N : cyclic_i;
}

inline int LETKF_worker::cycle(int i, int N)
{
	int cyclic_i = i % N;
	return cyclic_i<0 ? cyclic_i + N : cyclic_i;
}
		
inline void LETKF_worker::get_local(int middle, const vector &global, vector &local)
{
	for(int i=0; i<local.size(); ++i)
		local[i] = global[cycle(i + middle - span, global.size())];
}


// Local Ensemble Transform Kalman Filter for sub region

LETKF_region_master::LETKF_region_master
( int span, int span_avg, double epsilon, int region_begin, int region_end, int np, int rank ) 
: span(span), span_avg(span_avg), rho(1+epsilon), weights(0), region_begin(region_begin), region_end(region_end), np(np), rank(rank)
{
	set_weights();
}

LETKF_region_worker::LETKF_region_worker
( int span, int span_avg, double epsilon, int region_begin, int region_end, int np, int rank ) 
: span(span), span_avg(span_avg), rho(1+epsilon), weights(0), region_begin(region_begin), region_end(region_end), np(np), rank(rank)
{
	set_weights();
}

void LETKF_region_master::set_parameters( int span, int span_avg, double epsilon )
{
	LETKF_region_master::span = span;
	LETKF_region_master::span_avg = span_avg;
	rho = 1 + epsilon;
	set_weights();
}

void LETKF_region_worker::set_parameters( int span, int span_avg, double epsilon )
{
	LETKF_region_worker::span = span;
	LETKF_region_worker::span_avg = span_avg;
	rho = 1 + epsilon;
	set_weights();
}

void LETKF_region_master::set_weights()
{
	weights.resize( 2*span_avg + 1 );
	weights = 1. / weights.size();
}

void LETKF_region_worker::set_weights()
{
	weights.resize( 2*span_avg + 1 );
	weights = 1. / weights.size();
}

int LETKF_region_master::get_span()
{
	return span;
}

int LETKF_region_master::analysis_ensemble( matrix &Xg, const observation &obsv ) // for a sub region
{
	int flag = 0;
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
	
	//make Yg and compute xg_mean and yg_mean
	for( int i=0; i<k; ++i )
	{
		obsv.state_to_ob(Xg[i], Yg[i]);
		xg_mean += Xg[i];
		yg_mean += Yg[i];
	}	
	xg_mean /= k;
	yg_mean /= k;

	//subtract xg_mean and yg_mean from Xg and Yg
	for( int i=0; i<k; ++i ) {
		Xg[i] -= xg_mean;
		Yg[i] -= yg_mean;
	}

	for( int i=0; i<k; ++i ) {
		MPI::COMM_WORLD.Bcast( Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0 );
	}
	MPI::COMM_WORLD.Bcast( xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0 );

	// tag : locations ( 0 ~ n_data-1 ), end ( n_data )

	MPI::Status status;
	int buffer_ncols = 2 * span_avg + 1, buffer_size = k * buffer_ncols;
	double *buffer = new double[ buffer_size ];
	int n_data = region_end - region_begin + 1, n_sent = 0;
	int n_init_send = n_data > np-1 ? np-1 : n_data;
	
	// send min(n_data, np-1) data initially
	for( int i = 0; i < n_init_send; ++i, ++n_sent )
		MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, i+1, i );

	// receive and send more points 
	for( int sender, data_id, rec = 0; rec < n_data; ++rec ) // receive n_data times
	{
		MPI::COMM_WORLD.Recv( buffer, buffer_size, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
		sender = status.Get_source();
		
		if( n_sent < n_data ) {
			MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_sent++ );
		} else MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_data );

		data_id = status.Get_tag();
		if( data_id > n_data ) { flag = -1; continue; }
		
		// contribute to the global analysis ensemble		
		for( int count = 0, jg, j = 0; j < buffer_ncols; ++j )
		{
			jg = cycle( region_begin + data_id - span_avg + j, mg );
			for( int i = 0; i < k; ++i ) Xag[i][jg] += buffer[ count++ ];
		}
	}
	
	// adjust the boundary values
	double x;

	if( region_end - region_begin + 1 >= 2 * span_avg )
	{
		for( int jg1, jg2, j = 0; j < span_avg; ++j )
		{
			x = double( 2*span_avg + 1 ) / ( span_avg + 1 + j );
			jg1 = region_begin + j;
			jg2 = region_end - j;
			
			for( int i = 0; i < k; ++i ) { Xag[i][ jg1 ] *= x; Xag[i][ jg2 ] *= x; }
		}
	}
	else
	{
		for( int r, s, jg = region_begin; jg <= region_end; ++jg )
		{
			r = 1;
			s = jg - region_begin;
			r += s < span_avg ? s : span_avg;
			s = region_end - jg;
			r += s < span_avg ? s : span_avg;
			
			x = double( 2*span_avg + 1 ) / r;
			for( int i = 0; i < k; ++i ) Xag[i][ jg ] *= x;
		}
	}
		
	Xg = Xag;
	
	delete [] buffer;
	
	return flag;
}

void LETKF_region_worker::analysis_ensemble( matrix &Xg, const observation &obsv ) // for a sub region
{
	int flag;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	int n_data = region_end - region_begin + 1;
	
	matrix Yg(k, lg);
	vector xg_mean(mg), yg_mean(lg);
	
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
	vector temp( k );
	matrix Xa( k, m ); //local analysis ensemble
	
	for(int i=0; i<k; ++i) {
		MPI::COMM_WORLD.Bcast( Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0 );
	}
	MPI::COMM_WORLD.Bcast( xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0 );
	
	if( rank > n_data ) return;
	
	int buffer_ncols = 2 * span_avg + 1, buffer_size = k * buffer_ncols;
	double *buffer = new double[ buffer_size ];
	MPI::Status status;

	// compute local analysis ensemble for a local patch
	int x, data_id, n;
	while( 1 )
	{
		MPI::COMM_WORLD.Recv( &x, 0, MPI::INT, 0, MPI::ANY_TAG, status );
		data_id = status.Get_tag();
		if( data_id == n_data ) break;
		n = region_begin + data_id;
		
		// get local state and observation ensemble
		for( int i = 0; i < k; ++i )
		{
			get_local( n, Xg[i], X[i] );
			obsv.global_to_local( n, span, mg, Yg[i], Y[i] );
		}

		// get local x_mean, y_mean, yo
		get_local( n, xg_mean, x_mean );
		obsv.global_to_local( n, span, mg, yg_mean, y_mean );
		obsv.global_to_local( n, span, mg, obsv.get_observation(), yo );
		
		if( Y.size_col() == 0 ) { Xa = X; goto contribute; } 

		// compute C
		div( Y, obsv.covariance(), C );
		transpose( C );

		// compute Pw and W
		mul( Y, C, Pw );
		for( int i = 0; i < k; ++i ) Pw[i][i] += (k-1) / rho;
		flag = eigen( Pw, e, M );
		if( flag != 0 ) { data_id = n_data + 1; goto contribute; }
		flag = inverse( M, Minv );
		if( flag != 0 ) { data_id = n_data + 1; goto contribute; }
		for( int i=0; i < k; ++i ) for( int j = 0; j < k; ++j ) M[i][j] *= sqrt( (k-1) / e[j] );
		mul( M, Minv, W );
		flag = inverse( Pw );
		if( flag != 0 ) { data_id = n_data + 1; goto contribute; }
		
		// compute w_mean
		yo -= y_mean;
		mul( yo, C, temp );
		mul( temp, Pw, w_mean );

		// add w_mean to W
		for( int i=0; i < k; ++i ) W[i] += w_mean;

		// compute analysis ensemble
		mul( W, X, Xa );

		contribute:
		for( int i = 0; i < k; ++i ) Xa[i] += x_mean;
		
		double weight;
		for( int buffer_location = 0, jl, j = 0; j < buffer_ncols; ++j )
		{
			weight = weights[j];
			jl = span - span_avg + j;
			for( int i = 0; i < k; ++i ) buffer[ buffer_location++ ] = weight * Xa[i][jl];
		}
			
		// send the result
		MPI::COMM_WORLD.Send( buffer, buffer_size, MPI::DOUBLE, 0, data_id );
	}
	
	delete [] buffer;
}		

inline int LETKF_region_master::cycle( int i, int N )
{
	int cyclic_i = i % N;
	return cyclic_i < 0 ? cyclic_i + N : cyclic_i;
}

inline int LETKF_region_worker::cycle( int i, int N )
{
	int cyclic_i = i % N;
	return cyclic_i < 0 ? cyclic_i + N : cyclic_i;
}
		
inline void LETKF_region_worker::get_local( int middle, const vector &global, vector &local )
{
	for( int i = 0; i < local.size(); ++i )
		local[i] = global[ cycle( i + middle - span, global.size() ) ];
}


// Local Ensemble Transform Kalman Filter for global and sub regions with feedback strategy

LETKF_feedback_master::LETKF_feedback_master( int span_g, int span_avg_g, double epsilon_g, int span_r, int span_avg_r, double epsilon_r, int region_begin, int region_end, int np, int rank )
: span_g(span_g), span_avg_g(span_avg_g), rho_g(1+epsilon_g), span_r(span_r), span_avg_r(span_avg_r), rho_r(1+epsilon_r), weights_g(0), weights_r(0), region_begin(region_begin), region_end(region_end), np(np), rank(rank)
{
	set_weights();
}

LETKF_feedback_worker::LETKF_feedback_worker( int span_g, int span_avg_g, double epsilon_g, int span_r, int span_avg_r, double epsilon_r, int region_begin, int region_end, int np, int rank )
: span_g(span_g), span_avg_g(span_avg_g), rho_g(1+epsilon_g), span_r(span_r), span_avg_r(span_avg_r), rho_r(1+epsilon_r), weights_g(0), weights_r(0), region_begin(region_begin), region_end(region_end), np(np), rank(rank)
{
	set_weights();
}

void LETKF_feedback_master::set_weights()
{
	weights_g.resize( 2*span_avg_g + 1 );
	weights_r.resize( 2*span_avg_r + 1 );

	weights_g = 1. / weights_g.size();
	weights_r = 1. / weights_r.size();
}

void LETKF_feedback_worker::set_weights()
{
	weights_g.resize( 2*span_avg_g + 1 );
	weights_r.resize( 2*span_avg_r + 1 );

	weights_g = 1. / weights_g.size();
	weights_r = 1. / weights_r.size();
}

int LETKF_feedback_master::analysis_ensemble( matrix &Xg_g, matrix &Xg_r, const observation &obsv )
{
	int flag = 0;
	// dimensions
	int k = Xg_g.size_row(); // number of ensemble members
	int mg_g = Xg_g.size_col(), mg_r = Xg_r.size_col(); // dimensions of global states
	int lg = obsv.size(); // dimension of global observation
	int m_g = span_g * 2 + 1, m_r = span_r * 2 + 1; // dimensions of local states

	// for global computation
	matrix Xa_g( k, mg_g, 0. ), Xa_r( k, mg_r, 0. ); // global analysis ensembles
	matrix Yg_g( k ), Yg_r( k ); // global observation ensembles corresponding to global state ensembles
	vector xg_mean_g( mg_g, 0. ), xg_mean_r( mg_r, 0. ); // means of global state ensembles
	vector yg_mean_g( lg, 0. ), yg_mean_r( lg, 0. ); // means of global observation ensembles

	//make Yg and compute xg_mean and yg_mean
	for( int i=0; i<k; ++i )
	{
		obsv.state_to_ob( Xg_g[i], Yg_g[i] );
		obsv.state_to_ob( Xg_r[i], Yg_r[i] );
		xg_mean_g += Xg_g[i];
		xg_mean_r += Xg_r[i];
		yg_mean_g += Yg_g[i];
		yg_mean_r += Yg_r[i];
	}	
	xg_mean_g /= k;
	xg_mean_r /= k;
	yg_mean_g /= k;
	yg_mean_r /= k;

	//subtract xg_mean and yg_mean from Xg and Yg
	for( int i=0; i<k; ++i )
	{
		Xg_g[i] -= xg_mean_g;
		Xg_r[i] -= xg_mean_r;
		Yg_g[i] -= yg_mean_g;
		Yg_r[i] -= yg_mean_r;
	}

	for( int i=0; i<k; ++i )
	{
		MPI::COMM_WORLD.Bcast( Xg_g[i].get_v(), Xg_g[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Xg_r[i].get_v(), Xg_r[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg_g[i].get_v(), Yg_g[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg_r[i].get_v(), Yg_r[i].size(), MPI::DOUBLE, 0 );
	}
	MPI::COMM_WORLD.Bcast( xg_mean_g.get_v(), xg_mean_g.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( xg_mean_r.get_v(), xg_mean_r.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( yg_mean_g.get_v(), yg_mean_g.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( yg_mean_r.get_v(), yg_mean_r.size(), MPI::DOUBLE, 0 );

	// tag : locations in model true ( 0 ~ mg_r - 1 ), end ( mg_r ), error( mg_r+1 )

	MPI::Status status;
	int buffer_ncols_g = 2 * span_avg_g + 1;
	int buffer_ncols_r = 2 * span_avg_r + 1;
	int buffer_size = k * ( buffer_ncols_g + buffer_ncols_r );
	double *buffer = new double[ buffer_size ];
	int n_init_sent = 0;
	int q = mg_r / mg_g ;

	// send and receive data
	for( int receiver, n = 0; n < mg_r; ++n )
	{
		if( ! ( n % q == 0 || n >= region_begin && n <= region_end ) ) continue;
		
		if( n_init_sent == np - 1 )
		{
			MPI::COMM_WORLD.Recv( buffer, buffer_size, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
			save( status, buffer, buffer_ncols_g, buffer_ncols_r, k, mg_g, mg_r, Xa_g, Xa_r, flag );
			receiver = status.Get_source();
		}
		else receiver = ++n_init_sent;
		
		MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, receiver, n );
	}

	// receive the remaining results
	for( int n = 0; n < n_init_sent; ++n )
	{
		MPI::COMM_WORLD.Recv( buffer, buffer_size, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
		save( status, buffer, buffer_ncols_g, buffer_ncols_r, k, mg_g, mg_r, Xa_g, Xa_r, flag );
	}

	// send out a message signifying the end of data
	for( int i = 1; i < np; ++i ) MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, i, mg_r );
	
	// adjust the boundary values for regional analysis
	double x;

	if( region_end - region_begin + 1 >= 2 * span_avg_r )
	{
		for( int jg1, jg2, j = 0; j < span_avg_r; ++j )
		{
			x = double( 2*span_avg_r + 1 ) / ( span_avg_r + 1 + j );
			jg1 = region_begin + j;
			jg2 = region_end - j;
			
			for( int i = 0; i < k; ++i ) { Xa_r[i][ jg1 ] *= x; Xa_r[i][ jg2 ] *= x; }
		}
	}
	else
	{
		for( int r, s, jg = region_begin; jg <= region_end; ++jg )
		{
			r = 1;
			s = jg - region_begin;
			r += s < span_avg_r ? s : span_avg_r;
			s = region_end - jg;
			r += s < span_avg_r ? s : span_avg_r;
			
			x = double( 2*span_avg_r + 1 ) / r;
			for( int i = 0; i < k; ++i ) Xa_r[i][ jg ] *= x;
		}
	}
	
	Xg_g = Xa_g;
	Xg_r = Xa_r;
	
	delete [] buffer;

	return flag;
}

void LETKF_feedback_master::save( MPI::Status status, double *buffer, int buffer_ncols_g, int buffer_ncols_r, int k, int mg_g, int mg_r, matrix &Xa_g, matrix &Xa_r, int &flag )
{

	int id = status.Get_tag();
	if( id > mg_r ) { flag = -1; return; }
	int q = mg_r / mg_g, n_g = id / q, n_r = id, count = 0;

	if( id >= region_begin && id <= region_end && id % q == 0 ) // for global and regional
	{
		// contribute to the global analysis ensemble
		for( int jg, j = 0; j < buffer_ncols_g; ++j )
		{
			jg = cycle( n_g - span_avg_g + j, mg_g );
			for( int i = 0; i < k; ++i ) Xa_g[i][jg] += buffer[ count++ ];
		}
		
		// contribute to the regional analysis ensemble
		for( int jg, j = 0; j < buffer_ncols_r; ++j )
		{
			jg = cycle( n_r - span_avg_r + j, mg_r );
			for( int i = 0; i < k; ++i ) Xa_r[i][jg] += buffer[ count++ ];
		}
	}
	else if( id >= region_begin && id <= region_end && id % q != 0 ) // for regional
	{
		// contribute to the regional analysis ensemble
		for( int jg, j = 0; j < buffer_ncols_r; ++j )
		{
			jg = cycle( n_r - span_avg_r + j, mg_r );
			for( int i = 0; i < k; ++i ) Xa_r[i][jg] += buffer[ count++ ];
		}
	}
	else // for global
	{
		// contribute to the global analysis ensemble
		for( int jg, j = 0; j < buffer_ncols_g; ++j )
		{
			jg = cycle( n_g - span_avg_g + j, mg_g );
			for( int i = 0; i < k; ++i ) Xa_g[i][jg] += buffer[ count++ ];
		}
	}
}

void LETKF_feedback_worker::analysis_ensemble( matrix &Xg_g, matrix &Xg_r, const observation &obsv )
{
	int flag;
	// dimensions
	int k = Xg_g.size_row(); // number of ensemble members
	int mg_g = Xg_g.size_col(), mg_r = Xg_r.size_col(); // dimensions of global states
	int lg = obsv.size(); // dimension of global observation
	int m_g = span_g * 2 + 1, m_r = span_r * 2 + 1; // dimensions of local states
	double rho;
	
	matrix Yg_g( k, lg ), Yg_r( k, lg );
	vector xg_mean_g( mg_g ), xg_mean_r( mg_r ), yg_mean_g( lg ), yg_mean_r( lg );
		
	for( int i=0; i<k; ++i )
	{
		MPI::COMM_WORLD.Bcast( Xg_g[i].get_v(), Xg_g[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Xg_r[i].get_v(), Xg_r[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg_g[i].get_v(), Yg_g[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg_r[i].get_v(), Yg_r[i].size(), MPI::DOUBLE, 0 );
	}
	MPI::COMM_WORLD.Bcast( xg_mean_g.get_v(), xg_mean_g.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( xg_mean_r.get_v(), xg_mean_r.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( yg_mean_g.get_v(), yg_mean_g.size(), MPI::DOUBLE, 0 );
	MPI::COMM_WORLD.Bcast( yg_mean_r.get_v(), yg_mean_r.size(), MPI::DOUBLE, 0 );
	
	//for local computation
	matrix Xb_g( k, m_g ), Xb_r( k, m_r ); //local background ensemble
	matrix Xa_g( k, m_g ), Xa_r( k, m_r ); //local analysis ensemble
	matrix Y( k ); //local observation ensemble
	matrix C, Pw( k, k ), M( k, k ), Minv( k, k ), W( k, k );
	vector x_mean_g( m_g ), x_mean_r( m_r ), y_mean, yo, e( k ), w_mean( k ), temp( k );
	
	MPI::Status status;
	int buffer_ncols_g = 2 * span_avg_g + 1;
	int buffer_ncols_r = 2 * span_avg_r + 1;
	int buffer_size = k * ( buffer_ncols_g + buffer_ncols_r );
	double *buffer = new double[ buffer_size ];
	int dummy, id, n_g, n_r, q = mg_r / mg_g ;

	// compute local analysis ensemble for a local patch
	while( 1 )
	{
		MPI::COMM_WORLD.Recv( &dummy, 0, MPI::INT, 0, MPI::ANY_TAG, status );
		id = status.Get_tag();
		if( id == mg_r ) break;

		n_g = id / q, n_r = id;
		
		if( id >= region_begin && id <= region_end && id % q == 0 ) // for global and regional
		{
			for( int i = 0; i < k; ++i )
			{
				get_local( n_g, span_g, Xg_g[i], Xb_g[i] );
				get_local( n_r, span_r, Xg_r[i], Xb_r[i] );
				obsv.global_to_local( n_r, span_r, mg_r, Yg_r[i], Y[i] );
			}

			get_local( n_g, span_g, xg_mean_g, x_mean_g );
			get_local( n_r, span_r, xg_mean_r, x_mean_r );
			obsv.global_to_local( n_r, span_r, mg_r, yg_mean_r, y_mean );
			obsv.global_to_local( n_r, span_r, mg_r, obsv.get_observation(), yo );
			rho = rho_r;
		}
		else if( id >= region_begin && id <= region_end && id % q != 0 ) // for regional
		{
			for( int i = 0; i < k; ++i )
			{
				get_local( n_r, span_r, Xg_r[i], Xb_r[i] );
				obsv.global_to_local( n_r, span_r, mg_r, Yg_r[i], Y[i] );
			}

			get_local( n_r, span_r, xg_mean_r, x_mean_r );
			obsv.global_to_local( n_r, span_r, mg_r, yg_mean_r, y_mean );
			obsv.global_to_local( n_r, span_r, mg_r, obsv.get_observation(), yo );
			rho = rho_r;
		}
		else // for global
		{
			for( int i = 0; i < k; ++i )
			{
				get_local( n_g, span_g, Xg_g[i], Xb_g[i] );
				obsv.global_to_local( n_g, span_g, mg_g, Yg_g[i], Y[i] );
			}

			get_local( n_g, span_g, xg_mean_g, x_mean_g );
			obsv.global_to_local( n_g, span_g, mg_g, yg_mean_g, y_mean );
			obsv.global_to_local( n_g, span_g, mg_g, obsv.get_observation(), yo );
			rho = rho_g;
		}

		if( Y.size_col() == 0 ) { Xa_g = Xb_g; Xa_r = Xb_r; goto contribute; } 

		// compute C
		div( Y, obsv.covariance(), C );
		transpose( C );

		// compute Pw and W
		mul( Y, C, Pw );
		for( int i = 0; i < k; ++i ) Pw[i][i] += (k-1) / rho;
		flag = eigen( Pw, e, M );
		if( flag != 0 ) { id = mg_r+1; goto contribute; }
		flag = inverse( M, Minv );
		if( flag != 0 ) { id = mg_r+1; goto contribute; }
		for( int i=0; i < k; ++i ) for( int j = 0; j < k; ++j ) M[i][j] *= sqrt( (k-1) / e[j] );
		mul( M, Minv, W );
		flag = inverse( Pw );
		if( flag != 0 ) { id = mg_r+1; goto contribute; }
		
		// compute w_mean
		yo -= y_mean;
		mul( yo, C, temp );
		mul( temp, Pw, w_mean );

		// add w_mean to W
		for( int i=0; i < k; ++i ) W[i] += w_mean;

		// compute analysis ensemble
		if( id >= region_begin && id <= region_end && id % q == 0 ) // for global and regional
		{
			mul( W, Xb_g, Xa_g );
			mul( W, Xb_r, Xa_r );
		}
		else if( id >= region_begin && id <= region_end && id % q != 0 ) // for regional
			mul( W, Xb_r, Xa_r );
		else // for global
			mul( W, Xb_g, Xa_g );
		
		contribute:
		if( id >= region_begin && id <= region_end && id % q == 0 ) // for global and regional
		{
			for( int i = 0; i < k; ++i ) Xa_g[i] += x_mean_g;
			for( int i = 0; i < k; ++i ) Xa_r[i] += x_mean_r;
		}
		else if( id >= region_begin && id <= region_end && id % q != 0 ) // for regional
			for( int i = 0; i < k; ++i ) Xa_r[i] += x_mean_r;
		else // for global
			for( int i = 0; i < k; ++i ) Xa_g[i] += x_mean_g;

		double weight;
		int buffer_location = 0;

		if( id >= region_begin && id <= region_end && id % q == 0 ) // for global and regional
		{
			for( int jl, j = 0; j < buffer_ncols_g; ++j )
			{
				weight = weights_g[j];
				jl = span_g - span_avg_g + j;
				for( int i = 0; i < k; ++i ) buffer[ buffer_location++ ] = weight * Xa_g[i][jl];
			}
			
			for( int jl, j = 0; j < buffer_ncols_r; ++j )
			{
				weight = weights_r[j];
				jl = span_r - span_avg_r + j;
				for( int i = 0; i < k; ++i ) buffer[ buffer_location++ ] = weight * Xa_r[i][jl];
			}
			
			MPI::COMM_WORLD.Send( buffer, buffer_size, MPI::DOUBLE, 0, id );
		}
		else if( id >= region_begin && id <= region_end && id % q != 0 ) // for regional
		{
			for( int jl, j = 0; j < buffer_ncols_r; ++j )
			{
				weight = weights_r[j];
				jl = span_r - span_avg_r + j;
				for( int i = 0; i < k; ++i ) buffer[ buffer_location++ ] = weight * Xa_r[i][jl];
			}
			
			MPI::COMM_WORLD.Send( buffer, k * buffer_ncols_r, MPI::DOUBLE, 0, id );
		}
		else // for global
		{
			for( int jl, j = 0; j < buffer_ncols_g; ++j )
			{
				weight = weights_g[j];
				jl = span_g - span_avg_g + j;
				for( int i = 0; i < k; ++i ) buffer[ buffer_location++ ] = weight * Xa_g[i][jl];
			}
			
			MPI::COMM_WORLD.Send( buffer, k * buffer_ncols_g, MPI::DOUBLE, 0, id );
		}
	}

	delete [] buffer;
}		

inline int LETKF_feedback_master::cycle( int i, int N )
{
	int cyclic_i = i % N;
	return cyclic_i < 0 ? cyclic_i + N : cyclic_i;
}

inline int LETKF_feedback_worker::cycle( int i, int N )
{
	int cyclic_i = i % N;
	return cyclic_i < 0 ? cyclic_i + N : cyclic_i;
}
		
inline void LETKF_feedback_worker::get_local( int middle, int span, const vector &global, vector &local )
{
	for( int i = 0; i < local.size(); ++i )
		local[i] = global[ cycle( i + middle - span, global.size() ) ];
}


// Local Ensemble Transform Kalman Filter with bias
/*
LETKF_bias_master::LETKF_bias_master( int bias_number, int span, int span_avg, double epsilon, int np, int rank ) 
: bias_number(bias_number), span(span), span_avg(span_avg), rho(1+epsilon), weights(0), np(np), rank(rank)
{
	set_weights();
}

LETKF_bias_worker::LETKF_bias_worker( int span, int span_avg, double epsilon, int np, int rank ) 
: span(span), span_avg(span_avg), rho(1+epsilon), weights(0), np(np), rank(rank)
{
	set_weights();
}

void LETKF_bias_master::set_parameters( int span, int span_avg, double epsilon )
{
	LETKF_bias_master::span = span;
	LETKF_bias_master::span_avg = span_avg;
	rho = 1 + epsilon;
	set_weights();
}

void LETKF_bias_worker::set_parameters( int span, int span_avg, double epsilon )
{
	LETKF_bias_worker::span = span;
	LETKF_bias_worker::span_avg = span_avg;
	rho = 1 + epsilon;
	set_weights();
}

void LETKF_bias_master::set_weights()
{
	weights.resize( 2*span_avg + 1 );
	weights = 1. / weights.size();
}

void LETKF_bias_worker::set_weights()
{
	weights.resize( 2*span_avg + 1 );
	weights = 1. / weights.size();
}

int LETKF_bias_master::get_span()
{
	return span;
}

int LETKF_bias_master::analysis_ensemble( matrix &Xg, matrix &Bg, const observation &obsv )
{
	int flag = 0;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 +1; //dimension of local state
	MPI::Status status;
	
	//for global computation
	matrix Xag( k, mg, 0 ), Bag( k, mg, 0 ); //global analysis ensemble and global bias ensemble
	matrix Yg( k ); //global observation ensemble corresponding to global state ensemble
	vector xg_mean( mg, 0 ), bg_mean( mg, 0 ); //mean of global state ensemble and mean of global bias ensemble
	vector yg_mean( lg, 0 ); //mean of global observation ensemble
	
	//make Yg and compute xg_mean, bg_mean and yg_mean
	for( int i=0; i<k; ++i )
	{
		if( bias_number == 1 ) obsv.state_to_ob( Xg[i], Yg[i] );
		else if( bias_number == 2 ) obsv.state_to_ob( Xg[i] + Bg[i], Yg[i] );
		else { std::cerr << "unknown model bias number\n"; exit(-1); }
		xg_mean += Xg[i];
		bg_mean += Bg[i];
		yg_mean += Yg[i];
	}	
	xg_mean /= k;
	bg_mean /= k;
	yg_mean /= k;

	//subtract xg_mean, bg_mean and yg_mean from Xg, Bg and Yg
	for( int i=0; i<k; ++i )
	{
		Xg[i] -= xg_mean;
		Bg[i] -= bg_mean;
		Yg[i] -= yg_mean;
	}

	for( int i=0; i<k; ++i )
	{
		MPI::COMM_WORLD.Bcast( Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Bg[i].get_v(), Bg[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0 );
	}
	MPI::COMM_WORLD.Bcast( xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast( bg_mean.get_v(), bg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast( yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0);

	int n_data = mg / data_size;
	int n_init_send = n_data > np-1 ? np-1 : n_data;
	int n_sent = 0;
	
	// send min( n_data, np-1 ) data initially
	for( int i=0; i<n_init_send; ++i, ++n_sent)
		MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, i+1, i );

	int buffer_ncols = data_size + 2*span_avg;
	int buffer_size = 2 * k * buffer_ncols;
	double *buffer = new double[buffer_size];

	// receive and send more points 
	for( int sender, data_id, rec=0; rec<n_data; ++rec ) // receive n_data times
	{
		MPI::COMM_WORLD.Recv( buffer, buffer_size, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
		sender = status.Get_source();
		
		if( n_sent < n_data ) {
			MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_sent++ );
		} else MPI::COMM_WORLD.Send( MPI::BOTTOM, 0, MPI::INT, sender, n_data );

		data_id = status.Get_tag();
		if( data_id > n_data ) { flag = -1; continue; }
		
		//contribute to the global analysis ensemble		
		for( int count=0, jg, j=0; j<buffer_ncols; ++j )
		{
			jg = cycle( data_id * data_size - span_avg + j, mg );
			for( int i=0; i<k; ++i ) Xag[i][jg] += buffer[count++];
			for( int i=0; i<k; ++i ) Bag[i][jg] += buffer[count++];
		}
	}
		
	Xg = Xag;
	Bg = Bag;
	
	delete [] buffer;
	
	return flag;
}

void LETKF_bias_worker::analysis_ensemble( matrix &Xg, matrix &Bg, const observation &obsv)
{
	int flag;
	//dimensions
	int k = Xg.size_row(); //number of ensemble members
	int mg = Xg.size_col(); //dimension of global state
	int lg = obsv.size(); //dimension of global observation
	int m = span * 2 + 1; //dimension of local state
	int n_data = mg / data_size;
	
	matrix Yg( k, lg );
	vector xg_mean( mg ), bg_mean( mg ), yg_mean( lg );
	
	//for local computation
	matrix X( k, 2*m ); //local background ensemble
	matrix Y( k ); //local observation ensemble
	vector x_mean( 2*m );
	vector y_mean;
	vector yo; //local observation
	matrix C;
	matrix Pw( k, k );
	matrix M( k, k ); //for eigenvectors
	matrix Minv( k, k );
	vector e( k ); //for eigenvalues
	matrix W( k, k );
	vector w_mean( k );
	vector temp( k );
	matrix Xa( k, 2*m ); //local analysis ensemble

	for( int i=0; i<k; ++i )
	{
		MPI::COMM_WORLD.Bcast( Xg[i].get_v(), Xg[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Bg[i].get_v(), Bg[i].size(), MPI::DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( Yg[i].get_v(), Yg[i].size(), MPI::DOUBLE, 0 );
	}
	MPI::COMM_WORLD.Bcast( xg_mean.get_v(), xg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast( bg_mean.get_v(), bg_mean.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast( yg_mean.get_v(), yg_mean.size(), MPI::DOUBLE, 0);
	
	if ( rank > n_data ) return;
	
	int buffer_ncols = data_size + 2*span_avg;
	int buffer_size = 2 * k * buffer_ncols;
	double *buffer = new double[buffer_size];
	MPI::Status status;

	//compute local analysis ensemble for the partition
	int x, data_id, start_n, end_n;
	while(1)
	{
		MPI::COMM_WORLD.Recv( &x, 0, MPI::INT, 0, MPI::ANY_TAG, status );
		data_id = status.Get_tag();
		if( data_id==n_data ) break;
		
		for( int i=0; i<buffer_size; ++i ) buffer[i] = 0;
		int start_n = data_id * data_size;
		int end_n = start_n + data_size;
		for( int buffer_location_start=0, n=start_n; n<end_n; ++n, buffer_location_start += 2*k )
		{
			//get local state and observation ensemble
			for( int i=0; i<k; ++i )
			{
				get_local( n, Xg[i], Bg[i], X[i] );
				obsv.global_to_local( n, span, mg, Yg[i], Y[i] );
			}

			//get local x_mean, y_mean, yo
			get_local( n, xg_mean, bg_mean, x_mean );
			obsv.global_to_local( n, span, mg, yg_mean, y_mean );
			obsv.global_to_local( n, span, mg, obsv.get_observation(), yo );

			if( Y.size_col()==0 ) { Xa = X; goto contribute; } 

			//compute C
			div( Y, obsv.covariance(), C );
			transpose( C );

			//compute Pw and W
			mul( Y, C, Pw );
			for( int i=0; i<k; ++i ) Pw[i][i] += (k-1)/rho;
			flag = eigen( Pw, e, M );
			if( flag != 0 ) { data_id = n_data + 1; break; }
			flag = inverse( M, Minv );
			if( flag != 0 ) { data_id = n_data + 1; break; }
			for( int i=0; i<k; ++i ) for( int j=0; j<k; ++j ) M[i][j] *= sqrt((k-1)/e[j]);
			mul( M, Minv, W );
			flag = inverse( Pw );
			if( flag != 0 ) { data_id = n_data + 1; break; }
		
			//compute w_mean
			yo -= y_mean;
			mul( yo, C, temp );
			mul( temp, Pw, w_mean );

			//add w_mean to W
			for( int i=0; i<k; ++i ) W[i] += w_mean;

			//compute analysis ensemble
			mul( W, X, Xa );

			contribute:
			for( int i=0; i<k; ++i ) Xa[i] += x_mean;

			double weight;
			int avg_size = 2 * span_avg + 1;
			for( int buffer_location=buffer_location_start, jl, j=0; j<avg_size; ++j )
			{
				weight = weights[j];
				jl = span - span_avg + j;
				for( int i=0; i<k; ++i ) buffer[buffer_location++] += weight * Xa[i][jl];
				for( int i=0; i<k; ++i ) buffer[buffer_location++] += weight * Xa[i][jl+m];
			}
		}
		
		// send the result
		MPI::COMM_WORLD.Send( buffer, buffer_size, MPI::DOUBLE, 0, data_id );
	}
	
	delete [] buffer;
}		

inline int LETKF_bias_master::cycle(int i, int N)
{
	int cyclic_i = i % N;
	return cyclic_i<0 ? cyclic_i + N : cyclic_i;
}

inline int LETKF_bias_worker::cycle(int i, int N)
{
	int cyclic_i = i % N;
	return cyclic_i<0 ? cyclic_i + N : cyclic_i;
}
		
inline void LETKF_bias_worker::get_local( int middle, const vector &global1, const vector &global2, vector &local )
{
	int size_l = local.size()/2, size_g = global1.size();

	for( int i=0; i<size_l; ++i )
		local[i] = global1[ cycle( i + middle - span, size_g ) ];

	for( int i=0; i<size_l; ++i )
		local[size_l+i] = global2[ cycle( i + middle - span, size_g ) ];
}
*/

// Full Ensemble Kalman Filter

void FETKF_worker::analysis_ensemble(matrix &Xg, const observation &obsv) {}

