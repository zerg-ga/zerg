#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <vector>

class Derivative
{
public:
	Derivative();
	~Derivative();

	std::vector<double> Dfit(std::vector<double> &point, int type);

private:
	std::vector<double> DlennardJones(std::vector<double> &x);

	std::vector<double> Dgupta(std::vector<double> &x);


};

#endif
