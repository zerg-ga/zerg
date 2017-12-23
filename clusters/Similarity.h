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
		int method_in,
		int seed_in,
		int nAtoms_in,
		double tolSimilarity_in,
		double maxDistance_in,
		double minDistance_in,
		zerg::Printing * pPrinting_in,
		zerg::Random * rand_in);

	bool checkLimitations(std::vector<double> &x);
	bool checkSimilarity(std::vector<double> &x);
	bool checkSimilarity(std::vector<double> &x, std::vector< std::vector<double> > &targetIndividuals);
	void checkSimilarityGetRmsd(std::vector<double> &x);
	void appendTosimilarity();

	void addTargetIndividuals(
		std::vector< std::vector<double> > & x_vec,
		std::vector<int> & fitnessRank);

	void printNewBfgsInd();

private:
	double calcDistancesOverIAndGetMin(std::vector<double> &x, int i);
	std::vector<double> calcDistanceToCenter(std::vector<double> &x);
	std::vector<double> calcAndSortDistancesOverI(std::vector<double> &x, int i);
	std::vector<double> calcAndSortAllDistances(std::vector<double> &x);

	int seed;
	int method;
	int bestIndividualsSize;
	int printLevel; // 0-none ; 1-bfgs ; 2-all
	bool activateIntoBfgsRmsd;
	double minDistance;
	double maxDistance;
	double tolSimilarity;
	std::vector<double> tempDistance;
	std::vector<double> xTemp;
	std::vector< std::vector<double> > allDistances;
	std::vector< std::vector<double> > allCoordinates;
	std::vector< std::vector<double> > targetIndividuals;

	int nAtoms;

	zerg::Printing * pPrinting_;
	MarquesEnantiomers marqRmsd_;
	

};


#endif
