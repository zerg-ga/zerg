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
	int numberOfKilling;
	int maxGeneration;
	double highlanderInitialFitness;
	double highlanderThreshold;
	int highlanderMaxIteration;
	double adminLargeEnergyVariation;
	double adminMaxCreationVariation;
	int predatorMethod;
	double mutationValue;
	double crossoverWeight;
	double corssoverProbability;
	double gammaInitializeAtoms;
	double rcaInitializeAtoms;
	double maxDistance;
	double minDistance;
	int insistOnSimilar;
	int printingDebug;
	std::vector<double> initialCreationRate;

	// system parameters
	int numberOfParameters;
	int nAtomTypes1;
	int nAtomTypes2;
	int nAtomTypes3;
	int interactionPotentialType;
	int generationToChangeInteraction;
	std::vector<double> potentialParams;
	std::vector<int> atomTypes;
	std::vector<std::string> atomLabels;

	//operators parameters
	double scdo;
	double alfaMinGcdo;
	double alfaMaxGcdo;
	double wGcdo;
	double tetaMinTwisto;
	double tetaMaxTwisto;
	double contractionMinMtco;
	double contractionMaxMtco;

	//similarity parameters
	bool activateIntoBfgs;
	int similarityMethod;
	int similarityDebugLevel;
	double tolSimilarity;
	double energyReturnBfgs;
	int removeSimilarStructures;
	int bestIndividualSize;

	//util
	std::string experimentMethod;



};

}

#endif


