#include "Printing.h"

#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

using namespace std;

using namespace zerg;

namespace zerg {
	
	Printing::Printing()
	{
		mainOutput_.open("ga-output.txt");
	}

	Printing::~Printing()
	{
		mainOutput_.close();
	}

	void Printing::gaStartInputReading()
	{
		mainOutput_ << "Starting ga input reading..." << endl;
		mainOutput_ << "Lines:" << endl;
	}

	void Printing::showInputLines(string line)
	{
		mainOutput_ << line << endl;
	}

        void Printing::showAllParameters(
                zerg::GaParameters gaParam,
                string gamessNproc,
                string projectName,
                string interactionPotential,
                string gamessPath,
                string gamessScr,
                string gamessHeader,
                vector<string> baseFiles)
	{
		mainOutput_ << "END OF GA INPUT " << endl << endl << endl;
		mainOutput_ << "Printing all run information: " << endl;


		mainOutput_ << "Project name: " << projectName << endl;
		mainOutput_ << "Seed:  " << gaParam.seed << endl;
		mainOutput_ << "Restart: " << gaParam.restart << endl;
		mainOutput_ << "Number of cpus: " << gaParam.n_process << endl;
		mainOutput_ << "Population size: " << gaParam.pop_size << endl;
		mainOutput_ << "Max generation:  " << gaParam.maxGeneration << endl;
		mainOutput_ << "Highlander initial:  " << gaParam.highlanderInitialFitness << endl;
		mainOutput_ << "Highlander max generation: " << gaParam.highlanderMaxIteration << endl;
		mainOutput_ << "Admin energy variation: " << gaParam.adminLargeEnergyVariation << endl;
		mainOutput_ << "Admin max percentage: " << gaParam.adminMaxCreationVariation << endl;
		mainOutput_ << "Predator method: " << gaParam.predatorMethod << endl;
		mainOutput_ << "Mutation value: " << gaParam.mutationValue << endl;
		mainOutput_ << "CrossoverWeight: " << gaParam.crossoverWeight << endl;
		mainOutput_ << "Crossover probability: " << gaParam.corssoverProbability << endl;
		mainOutput_ << "Number of variables: " << gaParam.numberOfParameters << endl;
		mainOutput_ << "Gamma parameter: " << gaParam.gammaInitializeAtoms << endl;
		mainOutput_ << "rca parameter: " << gaParam.rcaInitializeAtoms << endl;
		mainOutput_ << "Max atom distance: " << gaParam.maxDistance << endl;
		mainOutput_ << "Min atom distance: " << gaParam.minDistance << endl;
		mainOutput_ << "Similar tries: " << gaParam.insistOnSimilar << endl;
		mainOutput_ << "Operators print debug: " << gaParam.printingDebug << endl;
		mainOutput_ << "scdo parameters: " << gaParam.scdo << endl;
		mainOutput_ << "alfaMinGcdo: " << gaParam.alfaMinGcdo << endl;
		mainOutput_ << "alfaMaxGcdo: " << gaParam.alfaMaxGcdo << endl;
		mainOutput_ << "wGcdo: " << gaParam.wGcdo << endl;
		mainOutput_ << "tetaMinTwisto: " << gaParam.tetaMinTwisto << endl;
		mainOutput_ << "tetaMaxTwisto: " << gaParam.tetaMaxTwisto << endl;
		mainOutput_ << "MinMtco: " << gaParam.contractionMinMtco << endl;
		mainOutput_ << "MaxMtco: " << gaParam.contractionMaxMtco << endl;
		for(size_t i = 0; i < gaParam.initialCreationRate.size(); i++)
		{
			mainOutput_ << "Creation rate " << i << ": " << gaParam.initialCreationRate[i] << endl;
		}

		mainOutput_ << "Interaction Potential: " << interactionPotential << endl;
		if(interactionPotential == "gamess")
		{
			mainOutput_ << "Gamess path: " << gamessPath << endl;
			mainOutput_ << "Gamess scr: " << gamessScr << endl;
			mainOutput_ << "Gamess header: " << gamessHeader << endl;
			mainOutput_ << "Number of bases: " << baseFiles.size() << endl;
			for(size_t j = 0; j < baseFiles.size(); j++)
			{
				mainOutput_ << "Base file " << j << ": " << baseFiles[j] << endl;
			}
		}
	}

	void Printing::endOfGamessOptions()
	{
		mainOutput_ << "Gamess options were sucessfully set" << endl;
	}

	void Printing::setCreationDebug(int creationDebug)
	{
		creationOutputDebugLevel = creationDebug;
		if (creationOutputDebugLevel == 1)
			creationOutput_.open("operators-histogram.txt");
		else if (creationOutputDebugLevel == 2)
			creationOutput_.open("operators-debug.txt");
	}	

	void Printing::showExperimentMethod(string experimentMethod)
	{
		mainOutput_ << "//////////////////////////////////////////" << endl
			<< "////////  ZERG GENETIC ALGORITHM  ////////" << endl
			<< "//////////////////////////////////////////" << endl
			<< endl << "Author: Frederico Teixeira Silva" << endl << endl;
		mainOutput_ << "Zerg GA started normally" << endl << endl;
		mainOutput_ << "Chosen GA option:  " << experimentMethod << endl << endl;
	}


	void Printing::writeOpenMessage()
	{
		mainOutput_ << "Input sucessfully read. Starting GA " << endl << endl;
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
			<< fixed << setprecision(16) << fitness_energies[0] << "   "
			<< fixed << setprecision(16) <<  fitness_energies[1] << "   "
			<< fixed << setprecision(16) <<  fitness_energies[2] << endl;

		mainOutput_ << "Dead individuals: " << endl
			<< dead_individuals[0] << "  "
			<< dead_individuals[1] << "  "
			<< dead_individuals[2] << "  "
			<< endl;

		mainOutput_ << "Fitness of the dead: " << endl
			<< fixed << setprecision(16) <<  fitness_energies[dead_individuals[0]] << "   "
			<< fixed << setprecision(16) <<  fitness_energies[dead_individuals[1]] << "   "
			<< fixed << setprecision(16) <<  fitness_energies[dead_individuals[2]] << endl
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
			mainOutput_ << endl << "highlander frequency:  " << frequency << endl << endl;
		}
	}

	void Printing::lastIndividualMessage(double frequency)
	{
		mainOutput_ << endl << "best individual frequency:  " << frequency << endl << endl;
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


