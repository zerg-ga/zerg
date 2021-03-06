#include "GeneticAlgorithm.h"

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <iomanip>

#include "Population.h"
#include "Predator.h"
#include "Creation.h"
#include "StructOptions.h"
#include "Printing.h"
#include "Random.h"

using namespace std;
using namespace zerg;

namespace zerg{
GeneticAlgorithm::GeneticAlgorithm(
	Population &pop_in,  
	GaParameters & gaParam,
	Random * rand_in,
	Printing * pPrinting_in)
:pop(pop_in)
{
	// tenho que seguir o geneticOut_
	// create histogram
	// era bom incluir o restart handling tambem
	pPrinting_ = pPrinting_in;
	pPrinting_->setCreationDebug(gaParam.printingDebug);
	pPrinting_->histogramTitle(gaParam.seed);

	generation = 1;
	maxGeneration = gaParam.maxGeneration;
	highlander = 0;
	pop_size = gaParam.pop_size;
	highlanderFitness = gaParam.highlanderInitialFitness;
	highlanderMaxIteration = gaParam.highlanderMaxIteration;
	highlanderThreshold = gaParam.highlanderThreshold;
	isToCheckMinimum = gaParam.isToCheckMinimum;

	//initializing objects
	pred_.initialize_predator(pop_size, gaParam.numberOfKilling, pPrinting_);
	creation_.initialize_creation(
		pop_size, 
		pop.get_number_of_creation_methods(),
		gaParam.n_process,
		pPrinting_,
		gaParam,
		rand_in);

	pPrinting_->writeOpenMessage();
}

GeneticAlgorithm::~GeneticAlgorithm()
{
	pPrinting_->endMessage();
}

void GeneticAlgorithm::ga_start()
{
	for (int i = 1; i <= maxGeneration; i++)
	{
		pPrinting_->generationMessage(i);
		predation();
		creation();

		if (checkHighlanderStop(i))
			return;
	}
	pPrinting_->generationEndMessage();
	double lastFreq = pop.checkMinimum(highlander);
	pPrinting_->lastIndividualMessage(lastFreq);
}

void GeneticAlgorithm::predation()
{
	pred_.get_dead_individuals(pop);
}

void GeneticAlgorithm::creation()
{
	creation_.make_new_individuals(pop,pred_);
}

bool GeneticAlgorithm::checkHighlanderStop(int i)
{
	for(int j=0; j<pop_size; j++)
	{
		if(pop.getfitness(j)<(highlanderFitness - highlanderThreshold))
		{
			highlander = j;
			highlanderFitness = pop.getfitness(j);
			highlanderFirstAppearence = i;
			pPrinting_->printHighlanderStatus(
					highlander,
					highlanderFirstAppearence,
					pop.getfitness(j));
		}
	}

	if((i - highlanderFirstAppearence) >= highlanderMaxIteration)
	{
		if(isToCheckMinimum){
			double firstFrequency = pop.checkMinimum(highlander);
			if (firstFrequency < 0.0e0)
			{
				pPrinting_->highlanderMessage(highlander, firstFrequency);
				highlanderFitness = 1.0e99;
				return false;	
			}
			pPrinting_->highlanderMessage(highlander, firstFrequency);
		}
		else
			pPrinting_->highlanderMessage(highlander, 0.0e0);
		
		return true;
	}
	else
		return false;
}

}
