#include "Printing.h"

#include <fstream>
#include <vector>

using namespace std;

using namespace zerg;

namespace zerg {
	
	Printing::Printing(int creationDebug)
	{
		mainOutput_.open("output-printing.txt");
		creationOutputDebugLevel = creationDebug;
		if (creationOutputDebugLevel == 1)
			creationOutput_.open("operators-histogram.txt");
		else if (creationOutputDebugLevel == 2)
			creationOutput_.open("operators-debug.txt");
	}

	Printing::~Printing()
	{
		mainOutput_.close();
	}

	void Printing::writeOpenMessage()
	{
		mainOutput_ << "//////////////////////////////////////////" << endl
			<< "////////  ZERG GENETIC ALGORITHM  ////////" << endl
			<< "//////////////////////////////////////////" << endl
			<< endl << "Author: Frederico Teixeira Silva" << endl << endl;

		mainOutput_ << "Zerg GA started normally" << endl << endl;

	}

	void Printing::generationMessage(int gen)
	{
		mainOutput_ << "Generation:  " << gen << endl;
	}

	void Printing::energyMessage(
		vector<int> &fitness_rank, 
		vector<int> &dead_individuals,
		vector<double> &fitness_energies)
	{
		mainOutput_ << "Best individuals: " << endl
			<< fitness_rank[0] << "   "
			<< fitness_rank[1] << "   "
			<< fitness_rank[2] << endl
			<< "Best fitness: " << endl
			<< fitness_energies[0] << "   "
			<< fitness_energies[1] << "   "
			<< fitness_energies[2] << endl;

		mainOutput_ << "Dead individuals: " << endl
			<< dead_individuals[0] << "  "
			<< dead_individuals[1] << "  "
			<< dead_individuals[2] << "  "
			<< endl;

		mainOutput_ << "Fitness of the dead: " << endl
			<< fitness_energies[dead_individuals[0]] << "   "
			<< fitness_energies[dead_individuals[1]] << "   "
			<< fitness_energies[dead_individuals[2]] << endl
			<< endl << endl;
	}

	void Printing::generationEndMessage()
	{
		mainOutput_ << "REACHED MAXIMUM GENERATIONS " << endl;
	}

	void Printing::highlanderMessage(int i, double frequency)
	{
		if (frequency < 0.0e0)
		{
			mainOutput_ << "Highlander "
				<< i
				<< " got a negative frequency (" << frequency << ")"
				<< " moving on."
				<< endl << endl;
		}
		else
		{
			mainOutput_ << "STOPPED BY HIGHLANDER SURVIVAL " << endl;
		}
	}

	void Printing::endMessage()
	{
		mainOutput_ << " Zerg GA terminated normally" << endl;
	}
	
	void Printing::similarityProblem(int method)
	{
		if (creationOutputDebugLevel == 2)
		{
			creationOutput_ << "WARNING!!!    Method:  "
				<< method << "   didnt surpass similarity. " << endl;
		}
	}

	void Printing::setNewIndividualsError()
	{
		if (creationOutputDebugLevel == 2)
		{
			creationOutput_ << " setNewIndividuals unexpected error " << endl;
		}
	}

	void Printing::variationOfEachMethod(int method, double variation)
	{
		if (creationOutputDebugLevel == 2)
		{
			creationOutput_ << "Method:  "
				<< method
				<< "   caused a variation on fitness of:  "
				<< variation
				<< endl;
		}
	}

	void Printing::allVariations(vector<double> &methodMean)
	{
		if (creationOutputDebugLevel == 2)
		{
			creationOutput_ << endl << "Mean fitness variations: " << endl;
			for (int ii = 0; ii < (int)methodMean.size(); ii++)
			{
				creationOutput_ << "Method: " << ii << "   mean:  " << methodMean[ii] << endl;
			}
			creationOutput_ << endl;
		}
	}

	void Printing::factorToIncreaseDecrease()
	{
		if(creationOutputDebugLevel == 2)
			creationOutput_ << "Factor to increase/decrease method N:  " << endl;
	}

	void Printing::printFactor(int method, double factor)
	{
		if (creationOutputDebugLevel == 2)
		{
			creationOutput_ << "method "
				<< method
				<< " : "
				<< factor
				<< endl;
		}
	}

	void Printing::normalizedCreationRate(vector<double> &creationRate)
	{
		if (creationOutputDebugLevel == 2)
		{
			creationOutput_ << endl << "Creation Rate after normalization " << endl;
			for (size_t i = 0; i < creationRate.size(); i++)
			{
				creationOutput_ << "method "
					<< i
					<< " creation rate:  "
					<< creationRate[i]
					<< endl;

			}
		}
	}

	void Printing::histogramTitle(int seed)
	{
		if(creationOutputDebugLevel > 0)
			creationOutput_ << "histogram of operators with seed:  " << seed << endl;
	}

	void Printing::histogramPrint(int method)
	{
		if (creationOutputDebugLevel > 0)
			creationOutput_ << method << "  ";
	}

	void Printing::histogramEndl()
	{
		if (creationOutputDebugLevel > 0)
			creationOutput_ << endl;
	}








}


