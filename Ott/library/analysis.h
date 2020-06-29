#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "tensor.h"
#include "model.h"
#include "observation.h"

// Ensemble transform Kalman filter

class ETKF
{
	char inflation_method;
	// 'r' for relaxation, 'm' for multiplicative in ensemble space, 'x' for multiplicative in state space
	double inflation;
	// inflation > 1 for a multiplicative covariance inflation
	// 0 < inflation < 1 weighting coefficient of background perturbations for relaxation method
public:
	ETKF( double epsilon, char inflation_method );
	int analysis( matrix &X, const observation &obsv );
};

class ETKF_regional
{
	int region_begin, region_end;
	char inflation_method;
	// 'r' for relaxation, 'm' for multiplicative in ensemble space, 'x' for multiplicative in state space
	double inflation;
	void get_region( const vector &global, vector &regional );
public:
	ETKF_regional( double epsilon, char inflation_method, int region_begin, int region_end );
	int analysis( matrix &Xg, const observation &obsv );
};

class ETKF_joint
{
	int region_begin, region_end;
	char inflation_method; // 'r' for relaxation, 'x' for multiplicative in state space
	double inflation, inflation_g, inflation_r;
	double Lambda; // y = Lambda * x_regional + (1-Lambda) * x_global
	double Kappa; // cost function = ... + Kappa * ( global - regional )^T * ( global - regional )
public:
	ETKF_joint( double epsilon_g, double epsilon_r, char inflation_method,
		int region_begin, int region_end, double Lambda, double Kappa );
	int analysis( matrix &Xg, matrix &Xr, const observation &obsv );
};

// Local ensemble transform Kalman filter

class LETKF
{
	int span; // number of grid points to the left and to the right from the local point
	vector weights; // for weighted average
	char inflation_method;
	// 'r' for relaxation, 'm' for multiplicative in ensemble space, 'x' for multiplicative in state space
	double inflation;
	// inflation > 1 for a multiplicative covariance inflation
	// 0 < inflation < 1 weighting coefficient of background perturbations for relaxation method
	void get_local( int middle, const vector &global, vector &local );
public:
	LETKF( int span, int span_avg, double epsilon, char inflation_method );
	LETKF( int span, const vector &weight_avg, double epsilon, char inflation_method );
	int analysis( matrix &Xg, const observation &obsv );
	int get_span();
};

class LETKF_regional
{
	int region_begin, region_end;
	int span;
	vector weights;
	char inflation_method;
	// 'r' for relaxation, 'm' for multiplicative in ensemble space, 'x' for multiplicative in state space
	double inflation;
	void get_local( int local_begin, int local_end, const vector &global, vector &local );
public:
	LETKF_regional( int span, int span_avg, double epsilon, char inflation_method, int region_begin, int region_end );
	int analysis( matrix &Xg, const observation &obsv );
};

class LETKF_joint
{
	int region_begin, region_end; // region_begin = 0, region_end = Nr-1 for the whole region without boundaries
	int span, span_avg; // in regional coordinate
	char inflation_method; // 'x' for multiplicative in state space
	double inflation_g, inflation_r;
	double Lambda; // y = Lambda * x_regional + (1-Lambda) * x_global
	double Kappa; // cost function = ... + Kappa * ( global - regional )^T * ( global - regional )
	void get_local( const vector &global, int begin, int end, vector &local, int start_local );
public:
	LETKF_joint( int span, int span_avg, double epsilon_g, double epsilon_r, char inflation_method, 
		int region_begin, int region_end, double Lambda, double Kappa );
	void analysis( matrix &Xg, matrix &Xr, const observation &obsv );
	//void analysis_cost( matrix &Xg, matrix &Xr, const observation &obsv, double &cost1, double &cost2, double &cost3 );
	//void bubble_sort( vector &v, matrix &m );
};

// Local ensemble Kalman filter

class LEKF
// based on 'a local ensemble Kalman filter for atmospheric data assimilation' Tellus 56A (2004), 5 E. Ott et al.
{
	int span; // number of grid points to the left and to the right from the local point
	int span_avg; // span for taking the average of analysis from neiboring points 
	double epsilon; // epsilon<1, covariance inflation factor
	void get_local( int middle, const vector &global, vector &local );
	void bubble_sort( vector &v, matrix &m );
public:
	LEKF( int span, int span_avg, double epsilon );
	int analysis( matrix &Xg_t, const observation &obsv, int k = 0 );
};

// Direct interpolation

class interpolation
{
public:
	void analysis( vector &x, const observation &obsv );
};

#endif
