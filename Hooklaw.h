#ifndef HOOKLAW_H
#define HOOKLAW_H

#include "BasicOperators.h"

namespace zerg{
class Hooklaw : public zerg::BasicOperators
{
public:
	Hooklaw(){};
	Hooklaw(int pop_size, int number_parameters);
	~Hooklaw();

private:
	void optimize(int ind_i);
	void start_hooklaw();

};

}

#endif

