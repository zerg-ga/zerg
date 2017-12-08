#include "Derivative.h"

#include <iostream>
#include <cmath>
#include <stdlib.h>

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

