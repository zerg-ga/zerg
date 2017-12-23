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
	FunctionDlibSim(int size_in, int fitType_in, Similarity * pSim_in)
	{
		size = size_in;
		fitType = fitType_in;
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

		return  fit_.fit(x, fitType);
	}

private:
	int size;
	int fitType;
	Similarity * pSim_;
};

#endif

#endif



