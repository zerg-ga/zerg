#include "BasicOperators.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include "Population.h"

using namespace std;
using namespace zerg;

namespace zerg{
BasicOperators::BasicOperators(int pop_size, int number_parameters)
{
	n_parameters = number_parameters;
	energy.resize(pop_size);
	x_vec.resize(pop_size);
	number_of_creation_methods = 5;
	tol_similarity = 1.0e-2;
	random_individual_range_max.resize(n_parameters);
	random_individual_range_min.resize(n_parameters);
}

BasicOperators::~BasicOperators(){}

void BasicOperators::startBasicOperators()
{
	mutationValue = 0.1e0;
	crossoverWeight = 0.7e0;
	crossoverProbability = 0.7e0;
}

void BasicOperators::startSimpleTest()
{
	int pop_size = energy.size();

	for(int j = 0; j<n_parameters; j++)
	{
		random_individual_range_min[j] = -10.e0;
		random_individual_range_max[j] = 10.0e0;
	}

	for(int i=0; i<pop_size; i++)
	{
		x_vec[i] = generate_random_individual();
		local_optimization(i);
	}
}

vector<double> BasicOperators::generate_random_individual()
{
	vector<double> x_rand(n_parameters);
	for (int i = 0; i<n_parameters; i++)
	{
		x_rand[i] = rand_->randomNumber(random_individual_range_min[i],
								random_individual_range_max[i]);
	}
	return x_rand;
}

void BasicOperators::optimize(int ind_i)
{
	int size = energy.size();
	double auxsoma = 0.0e0;
	for(int j=0; j<n_parameters; j++)
	{
		auxsoma += x_vec[ind_i][j]*x_vec[ind_i][j];
	}
	energy[ind_i] = auxsoma;
}


void BasicOperators::make_mutation(int target, int parent)
{
	for(int i=0; i<n_parameters; i++)
	{
		x_vec[target][i] = x_vec[parent][i];
		x_vec[target][i] += rand_->randomNumber(
			x_vec[target][i]-(x_vec[target][i]*mutationValue),
			x_vec[target][i]+(x_vec[target][i]*mutationValue));
	}
}

void BasicOperators::make_crossover_mean(int target, int parent1, int parent2)
{
	for(int i=0; i<n_parameters; i++)
	{
		x_vec[target][i] = crossoverWeight*(x_vec[parent1][i]+x_vec[parent2][i]);
	}
}

void BasicOperators::make_crossover_2_points(int target, int parent1, int parent2)
{
	int point1=0;
	int point2=0;
	select_2_points(point1,point2);

	for(int i=0;i<n_parameters;i++)
	{
		if((i>=point1)&&(i<=point2))
		{
			x_vec[target][i]=x_vec[parent2][i];
		}
		else
		{
			x_vec[target][i]=x_vec[parent1][i];
		}
	}
}

void BasicOperators::select_2_points(int &point1, int &point2)
{
	point1=-1;
	point2=-1;
	while(!(point2>point1))
	{
		point1 = rand_->randomNumber(0,n_parameters-2);
		point2 = rand_->randomNumber(point1+1,n_parameters-1);
	}
}

void BasicOperators::make_crossover_probability(int target,int parent1, int parent2)
{
	int size = (int) x_vec[target].size();
	double dice;
	for(int i=0;i<size;i++)
	{
		dice = rand_->randomNumber(0.0e0,1.0e0);
		if(dice>crossoverProbability)
			x_vec[target][i] = x_vec[parent1][i];
		else
			x_vec[target][i] = x_vec[parent2][i];
	}
}



////////
// GA //
////////

double BasicOperators::getfitness(int ind_i)
{
	return energy[ind_i];
}

void BasicOperators::local_optimization(int ind_i)
{
	optimize(ind_i);
	append_to_similarity(ind_i);
}

bool BasicOperators::create_individual(int creation_type,int target, int parent1, int parent2)
{
	switch(creation_type)
	{
	case 0:
		x_vec[target] = generate_random_individual();
		break;

	case 1:
		make_crossover_mean(target,parent1,parent2);
		break;

	case 2:
		make_crossover_2_points(target, parent1, parent2);
		break;

	case 3:
		make_mutation(target, parent1);
		break;

	case 4:
		make_crossover_probability(target, parent1, parent2);
		break;

	default:
		cout << "Creation type not avaible" << endl;
		return false;
	}
	return true;
}

bool BasicOperators::check_similarity(int target)
{
	return read_file_to_check_similarity(target);
}

void BasicOperators::append_to_similarity(int ind_i)
{
	write_similar_.open("similiarity.txt",fstream::app);
	write_similar_ << energy[ind_i] << endl;
	for(int i = 0; i<n_parameters; i++)
	{
		write_similar_ << x_vec[ind_i][i] << endl;
	}
	write_similar_.close();
}

bool BasicOperators::read_file_to_check_similarity(int target)
{
	read_similar_.open("similiarity.txt");
	double auxenergy;
	double auxx_veci;
	double auxerror;// least squares

	while(!read_similar_.eof())
	{
		read_similar_ >> auxenergy;
		auxerror = 0.0e0;
		for(int i=0;i<n_parameters;i++)
		{
			read_similar_ >> auxx_veci;
			auxx_veci-= x_vec[target][i];
			auxerror += abs(auxx_veci);
		}
		auxerror /= (double)n_parameters;

		if(auxerror<tol_similarity)
		{
			read_similar_.close();
			return true;
		}
	}
	read_similar_.close();
	return false;
}


bool BasicOperators::operatorAdministration(int method, const std::vector<double> &operatorPerformance)
{
	// if it has an administrator, pass to him.
	switch(method)
	{
	case 0:
		break;
	case 1:
		if(operatorPerformance[0] > 2.0e0)
			crossoverWeight = rand_->randomNumber(0.5e0,0.9e0);
		break;
	case 2:
		break;
	case 3:
		if(operatorPerformance[0] > 2.0e0)
			mutationValue = rand_->randomNumber(0.05e0,0.3e0);
		break;
	case 4:
		if(operatorPerformance[0] > 2.0e0)
			crossoverWeight = rand_->randomNumber(0.5e0,0.9e0);
		break;
	default:
		cout << " administration of this operator undefined contact developers " << endl;
		exit(2);
	}
	return true;
}


}
