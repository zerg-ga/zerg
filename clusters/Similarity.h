#ifndef SIMILARITY_H
#define SIMILARITY_H

#include <vector>
#include <string>

#include "InitializeAtoms.h"
#include "MarquesEnantiomers.h"
#include "../Printing.h"

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
		zerg::Printing * pPrinting_in);

	bool checkLimitations(std::vector<double> &x);
	bool checkSimilarity(std::vector<double> &x);
	bool checkSimilarity(std::vector<double> &x, std::vector< std::vector<double> > &targetIndividuals);
	std::vector<double>  checkSimilarityGetRmsd(std::vector<double> &x);
	void appendTosimilarity();

private:
	double calcDistancesOverIAndGetMin(std::vector<double> &x, int i);
	std::vector<double> calcDistanceToCenter(std::vector<double> &x);
	std::vector<double> calcAndSortDistancesOverI(std::vector<double> &x, int i);
	std::vector<double> calcAndSortAllDistances(std::vector<double> &x);

	int seed;
	int method;
	double minDistance;
	double maxDistance;
	double tolSimilarity;
	std::vector<double> tempDistance;
	std::vector< std::vector<double> > allDistances;

	int nAtoms;

	zerg::Printing * pPrinting_;
	MarquesEnantiomers marqRmsd_;

};


#endif
