#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include <vector>
#include <fstream>

#include "Population.h"
#include "Predator.h"
#include "Creation.h"
#include "StructOptions.h"
#include "Printing.h"

namespace zerg{	
class GeneticAlgorithm
{
public:
	GeneticAlgorithm(zerg::Population &Pop_in,GaParameters & gaParam);
	~GeneticAlgorithm();

	void ga_start();

private:
	int generation;
	int pop_size;
	int maxGeneration;
	
	bool checkHighlanderStop(int i);
	double highlanderFitness;
	int highlander;
	int highlanderMaxIteration;
	int highlanderFirstAppearence;

	zerg::Population &pop;
	zerg::Creation creation_;
	zerg::Predator pred_;

	void predation();
	void creation();

	zerg::Printing * pPrinting_;

};
}
#endif
