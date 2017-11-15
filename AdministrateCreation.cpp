#include "AdministrateCreation.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>

#include "FuzzyAdministration.h"

using namespace std;
using namespace zerg;

namespace zerg{
void AdministrateCreation::initializeAdministration(
	std::ofstream &geneticOut_, 
	GaOptions &gaoptions,
	GaParameters &gaParam)
{
	popSize = gaParam.pop_size;
	pgeneticOut_ = &geneticOut_;
	pgaoptions_ = &gaoptions;
	oldFitness.resize(popSize/4);
	methodUsed.resize(popSize/4);
	newIndividuals.resize(popSize/4);
	for(int i=0; i<popSize/4;i++)
		newIndividuals[i]=-1;

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
		*pgeneticOut_ << " setNewIndividuals unexpected error " << endl;
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

		if(pgaoptions_->printVariationOfCreationMethods)
		{
			*pgeneticOut_ << "Method:  " << methodUsed[i]
				<< "   caused a variation on fitness of:  " << varFitness
				<< endl;
		}
		newIndividuals[i] = -1;
	}
	if(pgaoptions_->printVariationOfCreationMethods)
	{
		*pgeneticOut_ << endl << "Mean fitness variations: " << endl;
		for(int ii=0;ii<(int)methodMean.size();ii++)
		{
			*pgeneticOut_ << "Method: " << ii << "   mean:  " << methodMean[ii] << endl;
		}
		*pgeneticOut_ << endl;
	}


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

	if(pgaoptions_->printVariationOfCreationMethods)
		*pgeneticOut_ << "Factor to increase/decrease method N:  " << endl;

	for(int i=0;i<nMethods;i++)
	{
		if(wasCreated(i,methodUsed))
		{
			creation_rate[i] += fuzzy_.getCreateRateVariation(methodMean[i]);
			if(pgaoptions_->printVariationOfCreationMethods)
			{
				*pgeneticOut_ << "method " 
					<< i 
					<< " : "
					<< fuzzy_.getCreateRateVariation(methodMean[i])
					<< endl;
			}
			if(creation_rate[i]<0.05)
				creation_rate[i]=0.05e0;
		}

	}

	//normalizing
	double auxsum = 0.0e0;
	for(int j=0;j<nMethods;j++)
	{
		auxsum+=creation_rate[j];
	}

	if(pgaoptions_->printVariationOfCreationMethods)
		*pgeneticOut_ << endl << "Creation Rate after normalization " << endl;
	for(int k=0;k<nMethods;k++)
	{
		creation_rate[k]/=auxsum;
		if(pgaoptions_->printVariationOfCreationMethods)
			*pgeneticOut_ 
			<< "method " << k << " creation rate:  " 
			<< creation_rate[k] 
			<< endl;
	}

	if(pgaoptions_->printVariationOfCreationMethods)
		*pgeneticOut_ << endl;
}

}

