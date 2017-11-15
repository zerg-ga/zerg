#ifndef STRUCTOPTIONS_H
#define STRUCTOPTIONS_H

#include <string>
#include <vector>

namespace zerg{
struct GaOptions
{
	bool printBestWorseIndividuals;
	bool printAllIndividuals;
	bool printVariationOfCreationMethods;
	bool similarityProblem;
};

// add type of fitness here
// n_process always 1
struct GaParameters
{
	int seed;
	bool restart;
	int n_process;
	int pop_size;
	int maxGeneration;
	double highlanderInitialFitness;
	int highlanderMaxIteration;
	double adminLargeEnergyVariation;
	double adminMaxCreationVariation;
	int predatorMethod;
	double mutationValue;
	double crossoverWeight;
	double corssoverProbability;
	int numberOfParameters;
	double gammaInitializeAtoms;
	double rcaInitializeAtoms;
	double maxDistance;
	int insistOnSimilar;
	std::vector<double> initialCreationRate;

	//operators parameters
	double scdo;
	double alfaMinGcdo;
	double alfaMaxGcdo;
	double wGcdo;
	double tetaMinTwisto;
	double tetaMaxTwisto;
	double contractionMinMtco;
	double contractionMaxMtco;

	//util
	std::string experimentMethod;



};

}

#endif


