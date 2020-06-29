#ifndef ANALYSIS_MPI_H
#define ANALYSIS_MPI_H

#include "mpi.h"
#include "tensor.h"
#include "model.h"
#include "observation.h"
#include "analysis.h"

// Local Ensemble Transform Kalman Filter

class LETKF_master {
	int span; // number of grid points to the left and to the right from the local point
	int span_avg; // span for taking the average of analysis from neiboring points 
	double rho; // rho>1, a multiplicative covariance inflation factor
	vector weights; // for weighted average
	int np, rank;
	void set_weights();
	int cycle(int i, int N);
public:
	LETKF_master(int span, int span_avg, double epsilon, int np, int rank);
	void set_parameters(int span, int span_avg, double epsilon);
	int get_span();
	int analysis_ensemble_old(matrix &Xg, const observation &obsv);
	int analysis_ensemble(matrix &Xg, const observation &obsv );
};

class LETKF_worker {
	int span; // number of grid points to the left and to the right from the local point
	int span_avg; // span for taking the average of analysis from neiboring points 
	double rho; // rho>1, a multiplicative covariance inflation factor
	vector weights; // for weighted average
	int np, rank;
	void set_weights();
	int cycle(int i, int N);
	void get_local(int middle, const vector &global, vector &local);
public:
	LETKF_worker(int span, int span_avg, double epsilon, int np, int rank);
	void set_parameters(int span, int span_avg, double epsilon);
	void analysis_ensemble_old(matrix &Xg, const observation &obsv);
	void analysis_ensemble(matrix &Xg, const observation &obsv );
};

// Local Ensemble Transform Kalman Filter for sub region
class LETKF_region_master {
	int span; // number of grid points to the left and to the right from the local point
	int span_avg; // span for taking the average of analysis from neiboring points 
	double rho; // rho>1, a multiplicative covariance inflation factor
	vector weights; // for weighted average
	int region_begin, region_end, np, rank;
	void set_weights();
	int cycle( int i, int N );
public:
	LETKF_region_master( int span, int span_avg, double epsilon, int region_begin, int region_end, int np, int rank );
	void set_parameters( int span, int span_avg, double epsilon );
	int get_span();
	int analysis_ensemble( matrix &Xg, const observation &obsv ); // for a sub region
};

class LETKF_region_worker {
	int span; // number of grid points to the left and to the right from the local point
	int span_avg; // span for taking the average of analysis from neiboring points 
	double rho; // rho>1, a multiplicative covariance inflation factor
	vector weights; // for weighted average
	int region_begin, region_end, np, rank;
	void set_weights();
	int cycle( int i, int N );
	void get_local( int middle, const vector &global, vector &local );
public:
	LETKF_region_worker( int span, int span_avg, double epsilon, int region_begin, int region_end, int np, int rank );
	void set_parameters( int span, int span_avg, double epsilon );
	void analysis_ensemble( matrix &Xg, const observation &obsv ); // for a sub region
};

// Local Ensemble Transform Kalman Filter for global and sub regions with feedback strategy
class LETKF_feedback_master {
	int span_g, span_r; // number of grid points to the left and to the right from the local point
	int span_avg_g, span_avg_r; // span for taking the average of analysis from neiboring points 
	double rho_g, rho_r; // rho>1, a multiplicative covariance inflation factor
	vector weights_g, weights_r; // for weighted average
	int region_begin, region_end, np, rank;
	void set_weights();
	int cycle( int i, int N );
	void save( MPI::Status status, double *buffer, int buffer_ncols_g, int buffer_ncols_r, int k, int mg_g, int mg_r, 
		matrix &Xa_g, matrix &Xa_r, int &flag );
public:
	LETKF_feedback_master( int span_g, int span_avg_g, double epsilon_g, int span_r, int span_avg_r, double epsilon_r, 
		int region_begin, int region_end, int np, int rank );
	int analysis_ensemble( matrix &Xg_g, matrix &Xg_r, const observation &obsv );
};

class LETKF_feedback_worker {
	int span_g, span_r; // number of grid points to the left and to the right from the local point
	int span_avg_g, span_avg_r; // span for taking the average of analysis from neiboring points 
	double rho_g, rho_r; // rho>1, a multiplicative covariance inflation factor
	vector weights_g, weights_r; // for weighted average
	int region_begin, region_end, np, rank;
	void set_weights();
	int cycle( int i, int N );
	void get_local( int middle, int span, const vector &global, vector &local );
public:
	LETKF_feedback_worker( int span_g, int span_avg_g, double epsilon_g, int span_r, int span_avg_r, double epsilon_r, 
		int region_begin, int region_end, int np, int rank );
	void analysis_ensemble( matrix &Xg_g, matrix &Xg_r, const observation &obsv );
};

// Local Ensemble Transform Kalman Filter with model bias
/*
class LETKF_bias_master {
	int bias_number; // model bias 1 or 2
	int span; // number of grid points to the left and to the right from the local point
	int span_avg; // span for taking the average of analysis from neiboring points 
	double rho; // rho>1, a multiplicative covariance inflation factor
	vector weights; // for weighted average
	int np, rank;
	void set_weights();
	int cycle( int i, int N );
public:
	LETKF_bias_master( int bias_number, int span, int span_avg, double epsilon, int np, int rank );
	void set_parameters( int span, int span_avg, double epsilon );
	int get_span();
	int analysis_ensemble( matrix &Xg, matrix &Bg, const observation &obsv );
};

class LETKF_bias_worker {
	int span; // number of grid points to the left and to the right from the local point
	int span_avg; // span for taking the average of analysis from neiboring points 
	double rho; // rho>1, a multiplicative covariance inflation factor
	vector weights; // for weighted average
	int np, rank;
	void set_weights();
	int cycle( int i, int N );
	void get_local( int middle, const vector &global1, const vector &global2, vector &local );
public:
	LETKF_bias_worker( int span, int span_avg, double epsilon, int np, int rank );
	void set_parameters( int span, int span_avg, double epsilon );
	void analysis_ensemble( matrix &Xg, matrix &Bg, const observation &obsv );
};
*/
// Ful Ensemble Transform Kalman Filter : dummy

typedef FETKF FETKF_master;

class FETKF_worker {
public:
	void analysis_ensemble(matrix &Xg, const observation &obsv);
};

#endif


