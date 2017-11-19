#include "GeneticAlgorithm.h"

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "Population.h"
#include "AuxMathGa.h"
#include "Predator.h"
#include "Creation.h"
#include "StructOptions.h"
#include "Printing.h"

using namespace std;
using namespace zerg;

// WARNING!!!
// pop_size has to be multiple of four.
namespace zerg{
GeneticAlgorithm::GeneticAlgorithm(
	Population &pop_in,  
	GaParameters & gaParam)
:pop(pop_in)
{
	// tenho que seguir o geneticOut_
	// create histogram
	// era bom incluir o restart handling tambem
	pPrinting_ = new Printing(1);
	pPrinting_->histogramTitle(gaParam.seed);

	generation = 1;
	maxGeneration = gaParam.maxGeneration;
	highlander = 0;
	pop_size = gaParam.pop_size;
	highlanderFitness = gaParam.highlanderInitialFitness;
	highlanderMaxIteration = gaParam.highlanderMaxIteration;
	geneticOut_.open("output.txt");

	//initializing objects
	pred_.initialize_predator(pop_size, pPrinting_, gaoptions_);
	creation_.initialize_creation(
		pop_size, 
		pop.get_number_of_creation_methods(),
		gaParam.n_process,
		geneticOut_,
		pPrinting_,
		gaoptions_,
		gaParam);

	setDefaultGaOptions();
	pPrinting_->writeOpenMessage();
	writeOpenMessage();
}

GeneticAlgorithm::~GeneticAlgorithm()
{
	pPrinting_->endMessage();
	delete pPrinting_;
	geneticOut_ << " Zerg GA terminated normally" << endl;
	geneticOut_.close();
}

void GeneticAlgorithm::ga_start()
{
	for (int i = 1; i <= maxGeneration; i++)
	{
		pPrinting_->generationMessage(i);
		geneticOut_ << "Generation:  " << i << endl;
		predation();
		creation();

		if (checkHighlanderStop(i))
			return;
	}
	pPrinting_->generationEndMessage();
}

void GeneticAlgorithm::setGaOptions(int flag, bool activate)
{
	switch (flag)
	{
	case 0:
		gaoptions_.printAllIndividuals = activate;
		break;

	case 1:
		gaoptions_.printBestWorseIndividuals = activate;
		break;

	case 2:
		gaoptions_.printVariationOfCreationMethods = activate;
		break;

	case 3:
		gaoptions_.similarityProblem = activate;
		break;

	default:
		cout << " option of GA not found " << endl;
		exit(2);
	}
}

void GeneticAlgorithm::writeOpenMessage()
{
	geneticOut_ << "//////////////////////////////////////////" << endl
				<< "////////  ZERG GENETIC ALGORITHM  ////////" << endl
				<< "//////////////////////////////////////////" << endl
				<< endl << "Author: Frederico Vultor" << endl << endl;

	geneticOut_ << "Zerg GA started normally" << endl << endl;
}

void GeneticAlgorithm::setDefaultGaOptions()
{
	setGaOptions(0,false);
	setGaOptions(1,true);
	setGaOptions(2,true);
	setGaOptions(3,true);
}

void GeneticAlgorithm::predation()
{
	pred_.get_dead_individuals(pop, geneticOut_);
}

void GeneticAlgorithm::creation()
{
	creation_.make_new_individuals(pop,pred_);
}

bool GeneticAlgorithm::checkHighlanderStop(int i)
{
	for(int j=0; j<pop_size; j++)
	{
		if(pop.getfitness(j)<highlanderFitness)
		{
			highlander = j;
			highlanderFitness = pop.getfitness(j);
			highlanderFirstAppearence = i;
		}
	}

	if((i - highlanderFirstAppearence)>highlanderMaxIteration)
	{
		pPrinting_->highlanderEndMessage();
		geneticOut_ << "STOPPED BY HIGHLANDER SURVIVAL " << endl;
		return true;
	}
	else
		return false;
}

}
