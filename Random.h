#ifndef RANDOM_H
#define RANDOM_H

#include <vector>

namespace zerg {
class Random
{
public:
	Random();

	~Random();

	double randomNumber(double fMin, double fMax);

	int randomNumber(int fMin, int fMax);

	std::vector<double> unitarySphericalVector();


private:


};
}

#endif