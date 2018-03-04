#include "Derivative.h"

#include <iostream>
#include <cmath>
#include <stdlib.h>

#include "Fitness.h"

using namespace std;

Derivative::Derivative(){}

Derivative::~Derivative(){}

vector<double> Derivative::Dfit(vector<double> &point, int type)
{
	switch (type)
	{
	case 0:
		return DlennardJones(point);
		break;
	case 1:
		cout << "ERROR ON: Fitness::Dfit - need to set parameters" << endl;
		exit(1);
		break;

	default:
		cout << "DERIVATIVA FUNCTION NOT FOUND" << endl;
		exit(3);
	}
}

vector<double> Derivative::Dfit(
	vector<double> &point, 
	int type,
	const vector<double> params,
	const vector<int> atomTypes)
{
	switch (type)
	{
	case 0:
		return DlennardJones(point);
		break;
	case 1:
		return Dgupta(point, params, atomTypes);
		break;

	default:
		cout << "DERIVATIVA FUNCTION NOT FOUND" << endl;
		exit(3);
	}
}


vector<double> Derivative::DlennardJones(vector<double> &x)
{
	// x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3
	int natm = x.size() / 3;
	double r, r2, r4, r8, r14;
	double dvlj;
	vector<double> dlj(3 * natm);
	for (int k = 0; k < 3 * natm; k++)
		dlj[k] = 0.0e0;

	for (int i = 0; i < (natm - 1); i++)
	{
		for (int j = (i + 1); j < natm; j++)
		{
			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
				);
			r2 = r * r;
			r4 = r2 * r2;
			r8 = r4 * r4;
			r14 = r8 * r4 * r2;

			dvlj = 4.0e0 * (6.0e0 / r8 - 12.0e0 / r14);
			dlj[i] += dvlj * (x[i] - x[j]);
			dlj[i + natm] += dvlj * (x[i + natm] - x[j + natm]);
			dlj[i + 2 * natm] += dvlj * (x[i + 2 * natm] - x[j + 2 * natm]);

			dlj[j] -= dvlj * (x[i] - x[j]);
			dlj[j + natm] -= dvlj * (x[i + natm] - x[j + natm]);
			dlj[j + 2 * natm] -= dvlj * (x[i + 2 * natm] - x[j + 2 * natm]);			
		}
	}
	return dlj;
}


vector<double> Derivative::Dgupta(
	vector<double> &x,
	const vector<double> params,
	const vector<int> atomTypes)
{
	vector<double> der(x.size());
	Fitness fit_;
	double e0 = fit_.fit(x, 1, params, atomTypes);
	double h = 1.0e-5;
	for (size_t i = 0; i < x.size(); i++)
	{
		vector<double> xh = x;
		xh[i] += h;
		double ei = fit_.fit(xh, 1, params, atomTypes);
		der[i] = (ei - e0) / h;
	}
	return der;
}




