#define useDlib

#ifdef useDlib

#ifndef FUNCTIONDLIBSIM_H
#define FUNCTIONDLIBSIM_H

#include <dlib/optimization.h>
#include <vector>

#include "Fitness.h"
#include "Similarity.h"

typedef dlib::matrix<double, 0, 1> column_vector;

class FunctionDlibSim
{
public:
	FunctionDlibSim(
		int size_in, 
		int fitType_in, 
		std::vector<double> &parameters,
		std::vector<int> &atomTypes_in,
		Similarity * pSim_in)
	{
		size = size_in;
		fitType = fitType_in;
		paramFunc = parameters;
		atomTypes= atomTypes_in;
		pSim_ = pSim_in;
	}

	~FunctionDlibSim(){}

	double operator() (const column_vector& arg) const
	{
		if (arg.size() == 0)
		{	
			return pSim_->checkSimilarityIntoBfgs();
		}

		std::vector<double> x(size);
		for (int i = 0; i < size; i++)
		{
			x[i] = arg(i);
		}
	
		pSim_->saveXToCheckBfgs(x);

		Fitness fit_;

		return  fit_.fit(x, fitType, paramFunc, atomTypes);
	}

private:
	int size;
	int fitType;
	std::vector<double> paramFunc;
	std::vector<int> atomTypes;
	Similarity * pSim_;
};

#endif

#endif



