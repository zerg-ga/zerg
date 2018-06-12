#ifndef SIMILARITY_H
#define SIMILARITY_H

#include <vector>
#include <string>

#include "InitializeAtoms.h"
#include "MarquesEnantiomers.h"
#include "Coordstructs.h"
#include "../Printing.h"
#include "../Random.h"

class Similarity
{
public:
	Similarity();
	~Similarity();

	void startSimilarity(
		zerg::GaParameters & gaParam,
		zerg::Printing * pPrinting_in,
		zerg::Random * rand_in);

	bool checkLimitations(std::vector<double> &x);
	int checkSimilarity(std::vector<double> &x);
	int checkSimilarity(std::vector<double> &x, std::vector< std::vector<double> > &targetIndividuals);
	void appendTosimilarity();

	double checkSimilarityIntoBfgs();
	void saveXToCheckBfgs(std::vector<double> &x);

	void addTargetIndividuals(
		std::vector< std::vector<double> > & x_vec,
		std::vector<int> & fitnessRank);

	void bestIndividualsCorrections(
        	std::vector<int> & corrections);

	void compareAllIndividuals(
        	std::vector< std::vector<double> > & x_vec,
        	std::vector<int> & fitnessRank);


	void printNewBfgsInd();

	void printBfgsSteps();

	void printEndIntoBfgs();
	
	void printInitialTitle();

	void printTitle();

	int printLevel; // 0-none ; 1-bfgs ; 2-all

private:
	double calcDistancesOverIAndGetMin(std::vector<double> &x, int i);
	std::vector<double> calcDistanceToCenter(std::vector<double> &x);
	std::vector<double> calcAndSortDistancesOverI(std::vector<double> &x, int i);
	std::vector<double> calcAndSortAllDistances(std::vector<double> &x);

	int nAtoms;
	int method;
	int bestIndividualsSize;
	int bfgsSteps;
	bool activateIntoBfgsRmsd;
	double minDistance;
	double maxDistance;
	double tolSimilarity;
	double energyReturnBfgs;
	std::vector<double> tempDistance;
	std::vector<double> xTemp;
	std::vector<double> xToCheckBfgs;
	std::vector< std::vector<double> > allDistances;
	std::vector< std::vector<double> > allCoordinates;
	std::vector< std::vector<double> > targetIndividuals;
	std::vector<std::string> atomLabels; 
	
	zerg::Printing * pPrinting_;
	MarquesEnantiomers marqRmsd_;
	

};


#endif
