
#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>

#include "RunAtomistica.h"

using namespace std;

int main()
{
	vector<double> x;
	x.push_back(0);
	x.push_back(1);
	x.push_back(2);
	x.push_back(0);
	x.push_back(0);
	x.push_back(0);
	x.push_back(0);
	x.push_back(0);
	x.push_back(0);
	double energy;
	RunAtomistica runAtom_;
	vector<double> optimized = runAtom_.run(x, energy);
	cout << "Expected energy:  -12.90" << endl;
	cout << "Calculated energy:  " << energy << endl;

}

