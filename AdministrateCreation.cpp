#include "AdministrateCreation.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>

#include "FuzzyAdministration.h"
#include "Printing.h"

using namespace std;
using namespace zerg;

namespace zerg{
void AdministrateCreation::initializeAdministration(
	zerg::Printing * pPrinting_in,
	GaParameters &gaParam)
{
	popSize = gaParam.pop_size;
	pPrinting_ = pPrinting_in;
	oldFitness.resize(popSize/4);
	methodUsed.resize(popSize/4);
	newIndividuals.resize(popSize/4);
	for(int i=0; i<popSize/4;i++)
		newIndividuals[i]=-1;

	fixed = gaParam.adminLargeEnergyVariation < 1.0e-5;
	fuzzy_.setFuzzyRules(gaParam.adminLargeEnergyVariation, gaParam.adminMaxCreationVariation);
}

void AdministrateCreation::setNewIndividuals(int newComer, int method, double fitness)
{
	int pos = -1;
	for(int i=0;i < (int)newIndividuals.size(); i++)
	{
		if(newIndividuals[i]==-1)
		{
			pos = i;
			break;
		}
	}

	if(pos==-1)
	{
		pPrinting_->setNewIndividualsError();
		exit(1);
	}

	newIndividuals[pos] = newComer;
	methodUsed[pos] = method;
	oldFitness[pos] = fitness;
}

void AdministrateCreation::adminCreationMethods(Population &pop, vector<double> &creation_rate)
{
	//setting variables
	// methodMean is mean of each method
	double varFitness;
	vector<double> methodMean(creation_rate.size());
	double nMethod = 1.0e0;
	int currentMethod = methodUsed[0];
	int size = (int)newIndividuals.size();
	for(int ii=0;ii<(int)creation_rate.size();ii++)
	{
		methodMean[ii] = 0.0e0;
	}

	double meanEnergy = 0.0e0;
	std::vector<int>::iterator it;
	for (int i = 0; i < popSize; i++)
	{
		it = find(newIndividuals.begin(), newIndividuals.end(), i);
		if( it == newIndividuals.end())
			meanEnergy += pop.getfitness(i);
	}
	meanEnergy /= (popSize - size);


	// calculating methodMean
	for(int i=0; i<size;i++)
	{
		varFitness = (pop.getfitness(newIndividuals[i]) - meanEnergy);
		methodMean[methodUsed[i]] += varFitness;
		if((currentMethod+1)!=methodUsed[i])
		{
			currentMethod+=1;
			methodMean[methodUsed[i]] /= nMethod;
			nMethod = 1.0e0;
		}
		else
		{
			nMethod+= 1.0e0;
		}

		pPrinting_->variationOfEachMethod(methodUsed[i], varFitness);
		newIndividuals[i] = -1;
	}
	pPrinting_->allVariations(methodMean);

	// administrate parameters of operators
	vector<double> operatorPerformance(1); // pass more info if needed
	for(int i=0; i<(int)methodMean.size();i++)
	{
		if(wasCreated(i,methodUsed))
		{
			operatorPerformance[0] = methodMean[i];
			pop.operatorAdministration(i,operatorPerformance);
		}
	}

	// modify
	setNewCreationRate(creation_rate,methodMean,methodUsed);
}

bool AdministrateCreation::wasCreated(int method, const vector<int> &methodUsed)
{
	for(int i=0; i<(int)methodUsed.size();i++)
	{
		if(methodUsed[i]==method)
			return true;
	}
	return false;
}

void AdministrateCreation::setNewCreationRate(vector<double> &creation_rate, const vector<double> &methodMean, const std::vector<int> &methodUsed)
{
	int nMethods = (int) methodMean.size();
	pPrinting_->factorToIncreaseDecrease();

	for(int i=0;i<nMethods;i++)
	{
		if(wasCreated(i,methodUsed))
		{
			creation_rate[i] += fuzzy_.getCreateRateVariation(methodMean[i]);
			pPrinting_->printFactor(i, fuzzy_.getCreateRateVariation(methodMean[i]));
			if((creation_rate[i]<0.05) && (!fixed))
				creation_rate[i]=0.05e0;
		}
	}

	//normalizing
	double auxsum = 0.0e0;
	for(int j=0;j<nMethods;j++)
	{
		auxsum+=creation_rate[j];
	}
	
	for(int k=0;k<nMethods;k++)
	{
		creation_rate[k]/=auxsum;
	}

	pPrinting_->normalizedCreationRate(creation_rate);
}

}

