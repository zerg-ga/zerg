#ifndef PREDATOR_H
#define PREDATOR_H

#include <vector>
#include <fstream>

#include "Population.h"
#include "StructOptions.h"
#include "Printing.h"

namespace zerg{
class Predator
{
public:
	void initialize_predator(
		int pop_size_in, 
		zerg::Printing * pPrinting_in);
	void get_dead_individuals(zerg::Population &pop);//change name to calculate

	std::vector<int> dead_individuals;

	std::vector<int> getFitnessRank();

	inline int get_pop_size(){return pop_size;}

private:
	int pop_size;
	std::vector<double> fitness_energies;
	std::vector<int> fitness_rank;

	zerg::Printing * pPrinting_;
};
}

#endif



