#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <string>

#include "Similarity.h"

class Fitness
{
public:
	Fitness();
	~Fitness();

	double fit(std::vector<double> &point, int type);

	double optimizeEmpiricalPotential(std::vector<double> &x, int fitType);

	double optimizeEmpiricalPotential(
		std::vector<double> &x, 
		int fitType,
		Similarity * pSim_);

	double runGamess(
		std::vector<double> &x,
		std::vector< std::string > &options,
		std::string gamessPath,
		std::string gamessScr,
		std::string nProc);

	double runGamessFrequency(
		int numberOfOptimizations,
		std::vector<double> &x,
		std::vector< std::string > &options,
		std::string gamessPath,
		std::string gamessScr,
		std::string nProc);


private:
	double lennardJones(std::vector<double> &x);

	double gupta(std::vector<double> &x);

	std::vector<double> getGuptaParameters(
		int i,
		int j,
		std::vector<int> & atomTypes,
		std::vector<double> & atomsParameters);

	double lowestPossibleEnergy;

};

#endif
