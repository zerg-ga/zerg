#include <iostream>
#include <string>

#include "GeneticAlgorithm.h"
#include "Printing.h"
#include "Random.h"

#include "clusters/ClustersFitness.h"
#include "clusters/ReadGaInput.h"

using namespace std;
using namespace zerg;

int main(int argc, char *argv[])
{
	string gaInput;
	if(argc == 3)
		gaInput = argv[2];
	else
		gaInput = "GaInput.txt";

	Printing printing_;
	ReadGaInput readGa_(&printing_);
	readGa_.inputName = gaInput;
	readGa_.readGaInput();
	srand(readGa_.getSeed());
	Random rand_;
	zerg::GaParameters gaParam = readGa_.getGaParameters();
	vector<string> options = readGa_.getOptions();
	ClustersFitness clFit_(
			&rand_,
			gaParam,
			options,
			readGa_.getGamessPath(),
			readGa_.getGamessScr(),
			readGa_.getGamessNprocess(),
			&printing_);
	GeneticAlgorithm ga1(
			clFit_, 
			gaParam,
			&rand_,
			&printing_);
	ga1.ga_start();
	clFit_.printAllIndividuals("finalPopulation.xyz");
	clFit_.printBfgsSteps();

	return 0;
}



