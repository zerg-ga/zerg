#include "Predator.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "AuxMathGa.h"
#include "Population.h"
#include "GAPrinting.h"

using namespace std;
using namespace zerg;

namespace zerg{
void Predator::initialize_predator(
	int pop_size_in,
	Printing * pPrinting_in)
{
	pop_size = pop_size_in;
	pPrinting_ = pPrinting_in;
	fitness_energies.resize(pop_size);
	fitness_rank.resize(pop_size);
	dead_individuals.resize(pop_size/4);
}

void Predator::get_dead_individuals(Population &pop)
{
	// ordering
	for(int i=0; i<pop_size; i++)
	{
		fitness_energies[i] = pop.getfitness(i);
		fitness_rank[i] = i;
	}
	vector<int> vector_order = AuxMathGa::vector_ordering(fitness_energies);
	AuxMathGa::vector_ordering_with_instructions(fitness_rank, vector_order);

	// choose last to kill
	int i_index = 0;
	for(int j=(pop_size-1); j>=((3*pop_size)/4); j--)
	{
		dead_individuals[i_index] = fitness_rank[j];
		i_index++;
	}

	pPrinting_->energyMessage(fitness_rank, dead_individuals,fitness_energies);


}

}