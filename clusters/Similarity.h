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
		bool activateIntoBfgsRmsd_in,
		int method_in,
		int seed_in,
		int nAtoms_in,
		int printLevel_in,
		double tolSimilarity_in,
		double maxDistance_in,
		double minDistance_in,
		zerg::Printing * pPrinting_in,
		zerg::Random * rand_in);

	bool checkLimitations(std::vector<double> &x);
	bool checkSimilarity(std::vector<double> &x);
	bool checkSimilarity(std::vector<double> &x, std::vector< std::vector<double> > &targetIndividuals);
	void appendTosimilarity();

	double checkSimilarityIntoBfgs();
	void saveXToCheckBfgs(std::vector<double> &x);

	void addTargetIndividuals(
		std::vector< std::vector<double> > & x_vec,
		std::vector<int> & fitnessRank);

	void printNewBfgsInd();

	void bestIndividualsCheck();

	void printBfgsSteps();

private:
	double calcDistancesOverIAndGetMin(std::vector<double> &x, int i);
	std::vector<double> calcDistanceToCenter(std::vector<double> &x);
	std::vector<double> calcAndSortDistancesOverI(std::vector<double> &x, int i);
	std::vector<double> calcAndSortAllDistances(std::vector<double> &x);

	int nAtoms;
	int seed;
	int method;
	int bestIndividualsSize;
	int printLevel; // 0-none ; 1-bfgs ; 2-all
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

	
	zerg::Printing * pPrinting_;
	MarquesEnantiomers marqRmsd_;
	

};


#endif
