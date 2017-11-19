#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <string>

class Fitness
{
public:
	Fitness();
	~Fitness();

	double fit(std::vector<double> &point, int type);

	double optimizeLennardJones(std::vector<double> &x, int fitType);

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

};

#endif
