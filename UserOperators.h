#ifndef USEROPERATORS_H
#define USEROPERATORS_H

#include <vector>

#include "BasicOperators.h"

class UserOperators : public zerg::BasicOperators
{
public:
	UserOperators(int pop_size, int number_parameters);
	~UserOperators();

	//reimplement or use basic
	void startUserOperators();//final adjustments  (like crossover probability)
	bool create_individual(int creation_type, int target, int parent1, int parent2);
	bool operatorAdministration(int method, const std::vector<double> &operatorPerformance);// modify operators with it's performance
	virtual bool check_similarity(int target) { return false; }
	virtual bool checkMinimum(int ind_i) { return true; }
	//keep
	virtual void local_optimization(int ind_i) = 0;
};

#endif



