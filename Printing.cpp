#include "Printing.h"

#define UNIX

#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <iostream>

using namespace std;

using namespace zerg;

namespace zerg {
	
Printing::Printing()
{
}

Printing::~Printing()
{
	mainOutput_.close();
	similarityOutput_.close();
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
	mainOutput_ << "Highlander threshold: " << gaParam.highlanderThreshold << endl;
	mainOutput_ << "Admin energy variation: " << gaParam.adminLargeEnergyVariation << endl;
	mainOutput_ << "Admin max percentage: " << gaParam.adminMaxCreationVariation << endl;
	mainOutput_ << "Predator method: " << gaParam.predatorMethod << endl;
	mainOutput_ << "Individuals predated: " << gaParam.numberOfKilling << endl;
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
	mainOutput_ << "activateIntoBfgs:  " << gaParam.activateIntoBfgs << endl;
	mainOutput_ << "Similarity method:  " << gaParam.similarityMethod << endl;
	mainOutput_ << "Similarity Debug level:  " << gaParam.similarityDebugLevel << endl;
	mainOutput_ << "Similirity tolerance:  " << gaParam.tolSimilarity << endl;
	mainOutput_ << "Energy return bfgs:  " << gaParam.energyReturnBfgs << endl;
	mainOutput_ << "Remove Similar Structures:  " << gaParam.removeSimilarStructures << endl;
	mainOutput_ << "Best individuals  size:  " << gaParam.bestIndividualSize << endl;
	for(size_t i = 0; i < gaParam.initialCreationRate.size(); i++)
	{
		mainOutput_ << "Creation rate " << i << ": " << gaParam.initialCreationRate[i] << endl;
	}
	mainOutput_ << "Interaction Potential: " << interactionPotential << endl;
	mainOutput_ << "Change interaction generation: " << gaParam.generationToChangeInteraction << endl;
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


void Printing::startingFirstIndividual()
{
	mainOutput_ << "Starting firt individual" << endl;
}

void Printing::endOfFirstIndividual()
{
	mainOutput_ << "First individual sucessful generated" << endl;
}

void Printing::endOfInitialPopulation()
{
	mainOutput_ << "Initial population sucessful generated" << endl;
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
	mainOutput_.open("ga-output.txt");
	mainOutput_ << "//////////////////////////////////////////" << endl
		<< "////////  ZERG GENETIC ALGORITHM  ////////" << endl
		<< "//////////////////////////////////////////" << endl
		<< endl << "Author: Frederico Teixeira Silva" << endl << endl;
	mainOutput_ << "Zerg GA started normally" << endl;

#ifdef UNIX
	time(&rawTime1);
	timeinfo = localtime(&rawTime1);
	char buffer[80];
	strftime(buffer, sizeof(buffer), "%d-%m-%Y %I:%M:%S",timeinfo);
	string dateI(buffer);
	mainOutput_ << "Date:  " << dateI << endl << endl;
#endif

	mainOutput_ << "Chosen GA option:  " << experimentMethod << endl << endl;
}

void Printing::openSimilarityFile()
{
	similarityOutput_.open("similarity-output.txt");
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
		<< scientific << dead_individuals[0] << "  "
		<< scientific << dead_individuals[1] << "  "
		<< scientific <<  dead_individuals[2] << "  "
		<< endl;

	mainOutput_ << "Fitness of the dead: " << endl
		<< fitness_energies[fitness_energies.size() - 1] << "   "
		<< fitness_energies[fitness_energies.size() - 2] << "   "
		<< fitness_energies[fitness_energies.size() - 3] << endl
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
	mainOutput_ << "Zerg GA terminated normally" << endl;
#ifdef UNIX
	time(&rawTime2);
	timeinfo = localtime(&rawTime2);
	char buffer[80];
	strftime(buffer, sizeof(buffer), "%d-%m-%Y %I:%M:%S",timeinfo);
	string dateI(buffer);
	mainOutput_ << "Ended at:  " << dateI << endl;
	mainOutput_ << "Wall time:  " << (int) rawTime2 - rawTime1
			<< " seconds" << endl;
#endif
}
	
void Printing::similarityProblem(int method)
{
	if (creationOutputDebugLevel == 2)
	{
		mainOutput_ << endl;
		mainOutput_ << "WARNING!!!    Method:  "
			<< method << "   didnt surpass similarity. " << endl;
		mainOutput_ << endl;
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

void Printing::printCreationIsEqual(double diff1, int i1)
{
	similarityOutput_ << "New individual is equal to "
		<< i1
		<< "  ->   value: "
		<< diff1
		<< endl
		<< endl;
}

void Printing::printSimilarityDistances(vector<double> &distances)
{
	for (size_t i = 0; i < distances.size(); i++)
	{
		similarityOutput_ << distances[i] << "  ";
	}
	similarityOutput_ << endl;
}

void Printing::printEqualStructures(vector<int> &equalStruct)
{
	similarityOutput_ << endl << "Equal structures:  ";
	for(size_t i = 0; i < equalStruct.size();i++)
		similarityOutput_ << equalStruct[i] << "  ";
	
	similarityOutput_ << endl;
}


void Printing::printSimilarityDistancesSelected(vector<double> &distances)
{
	similarityOutput_ << "selected structures" << endl;
	for (size_t i = 0; i < distances.size(); i++)
	{
		similarityOutput_ << distances[i] << "  ";
	}
	similarityOutput_ << endl;
}


void Printing::endlSimilarity()
{
	similarityOutput_ << endl 
			<< "STARTING INTO BFGS COMPARISSONS"
			<< endl << endl;
}

void Printing::printBfgsSteps(int bfgsSteps)
{
	mainOutput_ << endl
		<< "BFGS STEPS OF WHOLE SIMULATION"
		<< endl
		<< "steps: " << bfgsSteps
		<< endl << endl;
}

void Printing::printHighlanderStatus(
	int highI,
	int firstAppearence,
	double highFit)
{

	mainOutput_ << "BEST INDIVIDUAL CHANGE" << endl
		<< "best individual:  " << highI << endl
		<< "first appearance:  generation " << firstAppearence << endl
		<< "energy of best individual:  " << highFit << endl << endl;
}

void Printing::printSimInitOutBfgs()
{
	similarityOutput_ << "CALCULATING ROUTINE SIMILARITY" << endl;
}

void Printing::printSimUnique()
{
	similarityOutput_ << "INDIVIDUAL IS UNIQUE" << endl;
}

void Printing::printSimInitIntoBfgs()
{
	similarityOutput_ << "CALCULATING SIMILARITY INTO BFGS" << endl;
}

void Printing::printEndIntoBfgs()
{
	similarityOutput_ << endl
			<< "END OF INTO BFGS" << endl << endl;
}

void Printing::printSimTitle()
{
	similarityOutput_ << endl << endl
			<< "STARTING SIMILARITY" 
			<< endl << endl;
	
}

void Printing::printInitialSimTitle()
{
	similarityOutput_ << "FIRST POPULATION SIMILARITY" 
			<< endl << endl;
	
}



}


