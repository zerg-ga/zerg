#include "UserOperators.h"

#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace zerg;

UserOperators::UserOperators(int pop_size, int number_parameters)
:BasicOperators(pop_size, number_parameters)
{
//	number_of_creation_methods = 6;
//	tol_similarity = 1.0e-2;

}

UserOperators::~UserOperators(){}

void UserOperators::startUserOperators()
{
	startBasicOperators();
}

bool UserOperators::create_individual(int creation_type,int target, int parent1, int parent2)
{
	cout << "UserOperators::operatorAdministration is not active - contact developers" << endl;
	exit(1);
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

bool UserOperators::operatorAdministration(int method, const std::vector<double> &operatorPerformance)
{
	cout << "UserOperators::operatorAdministration is not active - contact developers" << endl;
	exit(1);
	// if it has an administrator, pass to him.
	switch(method)
	{
	case 0:
		break;
	case 1:
		//if(operatorPerformance[0] > 2.0e0)
			//crossoverWeight = AuxMathGa::randomNumber(0.5e0,0.9e0);
		break;
	case 2:
		break;
	case 3:
		//if(operatorPerformance[0] > 2.0e0)
			//mutationValue = AuxMathGa::randomNumber(0.05e0,0.3e0);
		break;
	case 4:
		//if(operatorPerformance[0] > 2.0e0)
			//crossoverWeight = AuxMathGa::randomNumber(0.5e0,0.9e0);
		break;
	default:
		cout << " administration of this operator undefined - contact developers " << endl;
		exit(2);
	}
	return true;
}





