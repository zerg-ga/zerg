#include "Hooklaw.h"

#include <iostream>
#include <vector>

#include "BasicOperators.h"

using namespace std;
using namespace zerg;

namespace zerg{
Hooklaw::Hooklaw(int pop_size, int number_parameters)
:BasicOperators(pop_size, number_parameters)
{
	start_hooklaw();
}

Hooklaw::~Hooklaw(){}

void Hooklaw::start_hooklaw()
{
	int pop_size = energy.size();

	// set creation range
	for(int j = 0; j<n_parameters; j++)
	{
		random_individual_range_min[j] = -1.e0;
		random_individual_range_max[j] = 1.0e0;
	}

	// can't touch
	bool aux;
	for(int i=0; i<pop_size; i++)
	{
		aux = create_individual(0,i,0,0); //method 0 always random
		local_optimization(i);
	}
	startBasicOperators();
}


void Hooklaw::optimize(int ind_i)
{
	// you have:
	// x_vec[ind_i] is a vector<double> -> to be optimized
	// I want:
	// energy[ind_i] -> fitness function 

	int size = energy.size();
	double auxsoma = 0.0e0;
	for(int j=0; j<n_parameters; j++)
	{
		auxsoma += x_vec[ind_i][j]*x_vec[ind_i][j];
	}
	energy[ind_i] = auxsoma;
}

}
