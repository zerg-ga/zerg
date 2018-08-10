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

	double fit(
		std::vector<double> &point, 
		int type,
		const std::vector<double> params,
		const std::vector<int> atomTypes);

	double optimizeEmpiricalPotential(
		std::vector<double> &x, 
		int fitType,
		std::vector<double> &parameters,
		std::vector<int> &atomTypes,
		Similarity * pSim_);

	double runGamess(
		std::vector<double> &x,
		std::vector< std::string > &options,
		std::string gamessPath,
		std::string gamessScr,
		std::string nProc);

	double runGamess(
		std::vector<double> &x,
		std::vector< std::string > &options,
		std::string gamessPath,
		std::string gamessScr,
		std::string nProc,
		Similarity * pSim_);

	double runGamessFrequency(
		int numberOfOptimizations,
		std::vector<double> &x,
		std::vector< std::string > &options,
		std::string gamessPath,
		std::string gamessScr,
		std::string nProc);

	double runMopac(
		std::vector<double> &x,
		std::vector<std::string> &options,
		std::vector<std::string> &atomLabels);


private:
	double lennardJones(std::vector<double> &x);

	double gupta(
		std::vector<double> &x,
		const std::vector<double> params,
		const std::vector<int> atomTypes);

	std::vector<double> getGuptaParameters(
		int i,
		int j,
		const std::vector<int> atomTypes,
		const std::vector<double> atomsParameters);

	double lowestPossibleEnergy;

	bool isFileExist(const std::string &name);

};

#endif
