#ifndef OBSERVATION_H
#define OBSERVATION_H

//#define SHIFT

#include <iostream>
#include <vector>
#include "tensor.h"
#include "math_lib.h"

class observation
{
	double sigma;	
	int n_ob;
	int *mark;
	// mark: contains indices of locations of observation in the true state at each grid point. 
	// It has -1 when there is no observation
	int *location;
	vector ob;
	int N; //size of a true state
	int cycle(int i, int N) const;
public:
	observation( double sigma, const std::vector<int> &locations, int size_true_state );
	~observation();
	int size_true_state() const;
	int *get_location();
	const int *get_location() const;
	int *get_mark() const;
	vector &get_observation();
	const vector &get_observation() const;
	double covariance() const;
	void cov_inv( int middle, int span, int size_ApprState, vector &r_inv ) const;
	void set_locations( const std::vector<int> &locations );
	void make_observation(const vector &TrueState);
	void state_to_ob( const vector &ApprState, vector &ob ) const;
	void Jacobian(const vector &ApprState, matrix &H) const;
	void ob_to_state(vector &ApprState) const;
	int size() const;
	void global_to_local
		(int middle, int span, int size_ApprState, const vector &ob_g, vector &ob_l) const;
	void global_to_local_Jacobian
		(int middle, int span, int size_ApprState, const matrix &Hg, matrix &H) const;
	void global_to_region( int region_begin, int region_end, const vector &ob_g, vector &ob_r ) const;
};

std::ostream& operator << (std::ostream& os, const observation &obsv);


////////////////////////////////// inline function definitions //////////////////////////////////////////////

inline int observation::cycle( int i, int N ) const { int j = i % N; return j<0 ? j + N : j; }

inline int observation::size() const { return n_ob; }

inline int observation::size_true_state() const { return N; }

inline int *observation::get_location() { return location; }

inline const int *observation::get_location() const { return location; }

inline int *observation::get_mark() const { return mark; }

inline vector &observation::get_observation() { return ob; }

inline const vector &observation::get_observation() const { return ob; }

inline double observation::covariance() const { return sigma * sigma; }

inline void observation::make_observation(const vector &TrueState)
{
	#ifdef SHIFT
	for( int i = n_ob/2; i < n_ob/2+n_ob; ++i )
		ob[ cycle(i,n_ob) ] = TrueState[ location[ cycle(i,n_ob) ] ] + random_gaussian( sigma );
	#else
	for( int i = 0; i < n_ob; ++i ) ob[i] = TrueState[ location[i] ] + random_gaussian( sigma );
	#endif
}

#endif

