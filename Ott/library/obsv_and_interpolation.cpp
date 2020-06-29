#include <algorithm>
#include "tensor.h"
#include "obsv_and_interpolation.h"

namespace
{
	inline int cycle(int i, int N)
	{
		int cyclic_i = i%N;
		return cyclic_i<0 ? cyclic_i+N : cyclic_i;
	}
}

void analysis_state(int n_ob, const vector &state, vector &a)
{
	int N = state.size();
	int selected_index;
	int temp;
	int *index = new int[N];

	for(int i=0; i<N; i++) index[i] = i;

	for(int i=0; i<n_ob; i++) {
        selected_index = i + int( (N-i) * ( rand() / (RAND_MAX+1.) ) );
		temp = index[i];
		index[i] = index[selected_index];
		index[selected_index] = temp;
	}

	std::sort(index, index+n_ob);

	int n0, n1, n2, n3;
	double z0, z1, z2, z3;
	int Lag0, Lag1, Lag2, Lag3;

	for(int i=0; i<n_ob; i++)
	{
		n0 = index[i];
		n1 = index[(i+1)%n_ob];
		n2 = index[(i+2)%n_ob];
		n3 = index[(i+3)%n_ob];

		z0 = state[n0];
		z1 = state[n1];
		z2 = state[n2];
		z3 = state[n3];

		n1 = n0 + cycle(n1-n0, N);
		n2 = n0 + cycle(n2-n0, N);
		n3 = n0 + cycle(n3-n0, N);

		a[n0] = state[n0];
		Lag0 = (n0-n1)*(n0-n2)*(n0-n3);
		Lag1 = (n1-n0)*(n1-n2)*(n1-n3);
		Lag2 = (n2-n0)*(n2-n1)*(n2-n3);
		Lag3 = (n3-n0)*(n3-n1)*(n3-n2);

		for(int n=n1+1; n<n2; n++) {
			a[cycle(n, N)] = static_cast<double>(z0*(n-n1)*(n-n2)*(n-n3))/Lag0
						+ static_cast<double>(z1*(n-n0)*(n-n2)*(n-n3))/Lag1
						+ static_cast<double>(z2*(n-n0)*(n-n1)*(n-n3))/Lag2
						+ static_cast<double>(z3*(n-n0)*(n-n1)*(n-n2))/Lag3;
		}
	}
	delete [] index;
}
