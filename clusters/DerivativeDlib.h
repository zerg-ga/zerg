#define useDlib

#ifdef useDlib

#ifndef DERIVATIVEDLIB_H
#define DERIVATIVEDLIB_H

#include <dlib/optimization.h>
#include <vector>

#include "Derivative.h"

typedef dlib::matrix<double, 0, 1> column_vector;

class DerivativeDlib
{
public:
	DerivativeDlib(
		int size_in, 
		int fitType_in,
		std::vector<double> & parameters,
		std::vector<int> & atomTypes_in)
	{
		size = size_in;
		fitType = fitType_in;
		paramFunc = parameters;
		atomTypes = atomTypes_in;
	}

	~DerivativeDlib(){}

	const column_vector operator() (const column_vector& arg) const
	{
		std::vector<double> x(size);
		for (int i = 0; i < size; i++)
		{
			x[i] = arg(i);
		}

		Derivative dfit_;

		std::vector<double> dxvector = dfit_.Dfit(
			x, 
			fitType, 
			paramFunc, 
			atomTypes);

		column_vector dxcolumn(size);
		for (int i = 0; i < size; i++)
		{
			dxcolumn(i) = dxvector[i];
		}

		return  dxcolumn;
	}

private:
	int size;
	int fitType;
	std::vector<double> paramFunc;
	std::vector<int> atomTypes;
};

#endif

#endif
