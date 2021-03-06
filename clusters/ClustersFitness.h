#ifndef CLUSTERSFITNESS_H
#define CLUSTERSFITNESS_H

#include <vector>
#include <string>
#include <fstream>

#include "ClustersOperators.h"
#include "../StructOptions.h"
#include "../Printing.h"
#include "../Random.h"

class ClustersFitness : public ClustersOperators
{
public:
	ClustersFitness(
		zerg::Random * rand_in,
		zerg::GaParameters & gaParam, 
		std::vector< std::string > &options_in,
		std::string gamessPath_in,
		std::string gamessScr_in,
		std::string nProc_in,
		zerg::Printing * pPrinting_in);
	~ClustersFitness();
	void local_optimization(int ind_i);
	void printAllIndividuals(std::string fileName);
	int getNumberOfLocalMinimizations() { return numberOfLocalMinimizations; }
	double checkMinimum(int ind_i);
	void printBfgsSteps();

private:
	void optimize(int ind_i);
	std::vector< std::string > options;
	std::vector< std::string > optionsGamess;
	std::string gamessPath;
	std::string gamessScr;
	std::string nProc;

	int numberOfLocalMinimizations;

	void saveIndividual(int ind_i);
	std::ofstream saveIndividual_;
	void readRestartFile();
	std::vector< std::vector<double> > restartCoordinates;
	std::vector<double> restartEnergies;
	int iRestart;
	int restartMax;

	// interation parameters
	int interactionType;
	std::vector<double> interactionParameters;
	std::vector<std::string> atomLabels;
	void changeInteraction(int newType);
	int generationToChangeInteraction;
	bool interactionChanged;

	// experiment util
	bool makeExperiment;
	double globalMinimaEnergy;
	double maxMinimizations;
	int seed;
	std::string experimentMethod;
	void setExperimentConditions(zerg::GaParameters & gaParam);
	void endExperimentConditions(int ind_i);


};

#endif


