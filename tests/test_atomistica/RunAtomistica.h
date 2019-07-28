#ifndef RUNATOMISTICA_H
#define RUNATOMISTICA_H

#include <vector>
#include <string>

class RunAtomistica
{
public:
	RunAtomistica(){};
	~RunAtomistica(){};

	std::vector<double> run(std::vector<double> & mol, double & energy);

private:
	void writeXyz(std::vector<double> & atoms, std::string testName);
	std::vector<double> readXyz(std::string xyzName, double & energy);

};

#endif

