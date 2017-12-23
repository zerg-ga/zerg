#include "Random.h"

#include <vector>
#include <cmath>
#include <stdlib.h>

using namespace std;
using namespace zerg;

Random::Random(){}

Random::~Random(){}

double Random::randomNumber(double fMin, double fMax)
{
	double f = ((double)rand() / (double)(RAND_MAX));
	return fMin + f * (fMax - fMin);
}

int Random::randomNumber(int fMin, int fMax)
{
	return fMin + (rand() % (int)(fMax - fMin + 1));
}

vector<double> Random::unitarySphericalVector()
{
	vector<double> unit(3);
	double pi = 3.1415926535897932384626433832795e0;
	double fi, teta;
	fi = 2.0e0 * pi * randomNumber(0.0e0, 1.0e0);
	teta = acos(2.0e0 * randomNumber(0.0e0, 1.0e0) - 1.0e0);
	unit[0] = sin(teta) * cos(fi);
	unit[1] = sin(teta) * sin(fi);
	unit[2] = cos(teta);
	return unit;
}




