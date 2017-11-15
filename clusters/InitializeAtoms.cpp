#include "InitializeAtoms.h"

#include <iostream>
#include <cmath>
#include <stdlib.h>

#include "../AuxMathGa.h"

using namespace std;

using namespace zerg;

InitializeAtoms::InitializeAtoms() 
{
        pi4sqr2 = 17.771531752633464988063523960243;
        pi = 3.1415926535897932384626433832795e0;
        exp3 = 0.33333333333333333333333333333333e0;
}

InitializeAtoms::~InitializeAtoms(){}

vector<double> InitializeAtoms::generateCluster(int natm_in, double gamma_in, double Rca_in)
{
	return generateClusterFred(natm_in, gamma_in, Rca_in);
}

vector<double> InitializeAtoms::generateClusterRondina(int natm_in, double gamma_in, double Rca_in)
{
	gamma = gamma_in;
	Rca = Rca_in;
	natm = natm_in;

	double rSphere = 2.0e0 * Rca *
		(0.5e0 + pow((3.0e0 * (double)natm) / pi4sqr2, exp3));

	vector<double> x(3 * natm);
	double r, teta, fi;
	double xi, yi, zi;
	int breakLoop = 0;
	do
	{
		for (int i = 0; i < natm; i++)
		{
			r = AuxMathGa::randomNumber(0.0e0, rSphere);
			//Sphere point picking to get equal area
			fi = 2.0e0 * pi * AuxMathGa::randomNumber(0.0e0, 1.0e0);
			teta = acos(2.0e0 * AuxMathGa::randomNumber(0.0e0, 1.0e0) - 1.0e0);
			xi = r * sin(teta) * cos(fi);
			yi = r * sin(teta) * sin(fi);
			zi = r * cos(teta);
			x[i] = xi;
			x[i + natm] = yi;
			x[i + 2 * natm] = zi;
		}
		breakLoop++;
		if (breakLoop > 1000000)
		{
			vector<double> xdummy;
			return xdummy;
		}
	} while (checkSuperposition(x, natm));

	return x;
}

vector<double> InitializeAtoms::generateClusterFred(int natm_in, double gamma_in, double Rca_in)
{
	gamma = gamma_in;
	Rca = Rca_in;
	natm = natm_in;

	double rSphere = 0.8e0 * Rca *
		(0.5e0 + pow((3.0e0 * (double)natm) / pi4sqr2, exp3));
	double lowerLimit = (1.0e0 - gamma) * Rca;
	vector<double> x(3 * natm);
	double r, teta, fi;
	double xi, yi, zi;
	int breakLoop = 0;
	for (int i = 0; i < natm; i++)
	{
		int k = 0;
		do
		{
			r = AuxMathGa::randomNumber(0.0e0, rSphere);
			//Sphere point picking to get equal area
			fi = 2.0e0 * pi * AuxMathGa::randomNumber(0.0e0, 1.0e0);
			teta = acos(2.0e0 * AuxMathGa::randomNumber(0.0e0, 1.0e0) - 1.0e0);
			xi = r * sin(teta) * cos(fi);
			yi = r * sin(teta) * sin(fi);
			zi = r * cos(teta);
			x[i] = xi;
			x[i + natm] = yi;
			x[i + 2 * natm] = zi;
			k++;
			if (k > 1000000)
			{
				cout << "infinity loop on Initialize atoms - check your parameters" << endl;
				exit(1);
			}
				
		} while (findMinimumDistanceUntilActualI(x, i) < lowerLimit);
	}

	return x;
}



double InitializeAtoms::calcDist(const vector<double> &x, int i, int j)
{
	return sqrt(
		(x[i] - x[j])*(x[i] - x[j]) +
		(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
		(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
		);
}

bool InitializeAtoms::checkSuperposition(const std::vector<double>& x, int iMax)
{
	double rMin;
	double lowerLimit = (1.0e0 - gamma) * Rca;
	double upperLimit = (1.0e0 + gamma) * Rca;
	for (int i = 0; i < iMax; i++)
	{
		rMin = findMinimumDistance(x, i);
		if ((rMin < lowerLimit) ||
			(rMin > upperLimit))
		{
			return true;
		}
	}

	return false;
}

double InitializeAtoms::findMinimumDistance(const std::vector<double>& x, int i)
{
	//distancia ao vizinho mais proximo
	double dist;
	double distmin = 1.0e3;
	for (int j = 0; j < natm; j++)
	{
		if (i != j)
		{
			dist = calcDist(x, i, j);
			if (dist < distmin)
			{
				distmin = dist;
			}
		}
	}
	return distmin;
}


double InitializeAtoms::findMinimumDistanceUntilActualI(const std::vector<double>& x, int i)
{
	//distancia ao vizinho mais proximo
	double dist;
	double distmin = 1.0e3;
	for (int j = 0; j < i; j++)
	{
		if (i != j)
		{
			dist = calcDist(x, i, j);
			if (dist < distmin)
			{
				distmin = dist;
			}
		}
	}
	return distmin;
}

