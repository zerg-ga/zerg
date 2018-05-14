#include "Random.h"

#include <vector>
#include <cmath>
#include <stdlib.h>

using namespace std;
using namespace zerg;

Random::Random()
{
	iff = 0;
}

Random::~Random(){}

double Random::randomNumber(double fMin, double fMax)
{
	double f = ((double)rand() / (double)(RAND_MAX));
	return fMin + f * (fMax - fMin);
}

double Random::randomNumber3(double fMin, double fMax)
{
	int idum = 1;
	double f = ((double)ran3(&idum));
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

vector<double> Random::unitaryCartesianVector()
{
	vector<double> unit(3);
	double x = randomNumber(-1.0e0,1.0e0);
	double y = randomNumber(-1.0e0,1.0e0);
	double z = randomNumber(-1.0e0,1.0e0);
	double r = sqrt(x*x + y*y + z*z);
	unit[0] = x / r;
	unit[1] = y / r;
	unit[2] = z / r;
	return unit;
}







#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float Random::ran3(int *idum)
{
	long mj, mk;
	int i, ii, k;

	if (*idum < 0 || iff == 0) {
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ) mk += MBIG;
			mj = ma[ii];
		}
		for (k = 1; k <= 4; k++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if (++inext == 56) inext = 1;
	if (++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return mj * FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC





