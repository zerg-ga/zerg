#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <vector>

class Derivative
{
public:
	Derivative();
	~Derivative();

	std::vector<double> Dfit(std::vector<double> &point, int type);

	std::vector<double> Dfit(
		std::vector<double> &point, 
		int type,
		const std::vector<double> params,
		const std::vector<int> atomTypes);

private:
	std::vector<double> DlennardJones(std::vector<double> &x);

	std::vector<double> Dgupta(
		std::vector<double> &x,
		const std::vector<double> params,
		const std::vector<int> atomTypes);


};

#endif
