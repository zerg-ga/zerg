#ifndef INITIALIZEATOMS_H
#define INITIALIZEATOMS_H

#include <vector>
#include <string>

#include "../Random.h"

class InitializeAtoms
{
public:
	InitializeAtoms();

	~InitializeAtoms();

	void initializeSetSeed(zerg::Random * rand_in);

	std::vector<double> generateCluster(int Natoms, double gamma_in, double Rca_in);
	
	double calcDist(const std::vector<double> &x, int i, int j);

private:
	double gamma;//parameter - sugerido: entre 0,2 e 0,4. Sao gerados entre (1-gamma) ate (1+gamma)
	double Rca;//Rcalfa na tese pag 62 - media do raio covalente.
	int natm;

	double pi4sqr2;

	double pi;

	double exp3;

	std::vector<double> generateClusterRondina(int Natoms, double gamma_in, double Rca_in);

	std::vector<double> generateClusterFred(int Natoms, double gamma_in, double Rca_in);

	bool checkSuperposition(const std::vector<double> &x, int iMax);

	double findMinimumDistance(const std::vector<double> &dist, int i);

	double findMinimumDistanceUntilActualI(const std::vector<double> &dist, int i);

	zerg::Random * rand_;
};

#endif



