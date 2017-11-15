#ifndef USERFITNESS_H
#define USERFITNESS_H

#include "UserOperators.h"

class UserFitness : public UserOperators
{
public:
	UserFitness(int pop_size, int number_parameters);

	~UserFitness();

	void local_optimization(int ind_i);

private:
	void optimize(int ind_i);
};

#endif

