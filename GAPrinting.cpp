#include "GAPrinting.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;


namespace zerg{
void GAPrinting::print_vector(const vector<int> &vec_in)
{
	cout << setiosflags(ios::fixed) << setw(12);

	int len = vec_in.size();
	for(int i = 0; i<len ; i++)
	{
		cout << vec_in[i] << "     ";
	}
	cout << endl;
}

void GAPrinting::print_vector(const vector<int> &vec_in, ofstream &currentOut)
{
	currentOut << setiosflags(ios::fixed) << setw(12);
	int len = vec_in.size();
	for(int i = 0; i<len ; i++)
	{
		currentOut << vec_in[i] << "     ";
	}
	currentOut << endl;
}


void GAPrinting::print_vector(const vector<double> &vec_in)
{
	cout << setiosflags(ios::fixed) <<  setw(12);

	int len = vec_in.size();
	for(int i = 0; i<len ; i++)
	{
		cout << vec_in[i] << "     ";
	}
	cout << endl;

}

void GAPrinting::print_vector(const vector<double> &vec_in, ofstream &currentOut)
{
	currentOut << setiosflags(ios::fixed) << setw(12);

	int len = vec_in.size();
	for(int i = 0; i<len ; i++)
	{
		currentOut << vec_in[i] << "     ";
	}
	currentOut << endl;
}


}