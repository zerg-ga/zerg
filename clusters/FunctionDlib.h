#define useDlib

#ifdef useDlib

#ifndef FUNCTIONDLIB_H
#define FUNCTIONDLIB_H

#include <dlib/optimization.h>
#include <vector>

#include "Fitness.h"

typedef dlib::matrix<double, 0, 1> column_vector;

class FunctionDlib
{
public:
	FunctionDlib(int size_in, int fitType_in)
	{
		size = size_in;
		fitType = fitType_in;
	}

	~FunctionDlib(){}

	double operator() (const column_vector& arg) const
	{
		std::vector<double> x(size);
		for (int i = 0; i < size; i++)
		{
			x[i] = arg(i);
		}

		Fitness fit_;

		return  fit_.fit(x, fitType);
	}

private:
	int size;
	int fitType;
};

#endif

#endif



