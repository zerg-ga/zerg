#ifndef PREDATOR_H
#define PREDATOR_H

#include <vector>
#include <fstream>

#include "Population.h"
#include "StructOptions.h"

namespace zerg{
class Predator
{
public:
	void initialize_predator(int pop_size_in, zerg::GaOptions &gaoptions);
	void get_dead_individuals(zerg::Population &pop, std::ofstream &geneticOut_);//change name to calculate

	std::vector<int> dead_individuals;
	inline int get_pop_size(){return pop_size;}

private:
	int pop_size;
	std::vector<double> fitness_energies;
	std::vector<int> fitness_rank;

	void printPredatorInfo(std::ofstream &geneticOut, 
                           const std::vector<int> &fitness_rank, 
						   const std::vector<double> &fitness_energies);

	zerg::GaOptions * pgaoptions_;
};
}

#endif



