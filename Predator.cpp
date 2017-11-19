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
	Printing * pPrinting_in,
	GaOptions &gaoptions)
{
	pgaoptions_ = &gaoptions;
	pop_size = pop_size_in;
	pPrinting_ = pPrinting_in;
	fitness_energies.resize(pop_size);
	fitness_rank.resize(pop_size);
	dead_individuals.resize(pop_size/4);
}

void Predator::get_dead_individuals(Population &pop, ofstream &geneticOut_)
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

	printPredatorInfo(geneticOut_, fitness_rank, fitness_energies);

	pPrinting_->energyMessage(fitness_rank, dead_individuals,fitness_energies);


}

// printing
void Predator::printPredatorInfo(ofstream &geneticOut_, 
								 const std::vector<int> &fitness_rank, 
								 const std::vector<double> &fitness_energies)
{
	if(pgaoptions_->printBestWorseIndividuals)
	{
		geneticOut_	<< "Best individuals: " << endl
			<< fitness_rank[0] << "   " 
			<<fitness_rank[1] << "   " 
			<<fitness_rank[2] << endl
			<< "Best fitness: " << endl
			<< fitness_energies[0] << "   "
			<<fitness_energies[1] << "   " 
			<<fitness_energies[2] << endl;

		geneticOut_ << "Dead individuals: " << endl
			<< dead_individuals[0] << "  "
			<< dead_individuals[1] << "  "
			<< dead_individuals[2] << "  "
			<< endl;

		geneticOut_ << "Fitness of the dead: " << endl
			<< fitness_energies[dead_individuals[0]] << "   "
			<< fitness_energies[dead_individuals[1]] << "   " 
			<< fitness_energies[dead_individuals[2]] << endl
			<< endl << endl;
	}

	if(pgaoptions_->printAllIndividuals)
	{
		geneticOut_ << "All individuals: " << endl;
		GAPrinting::print_vector(fitness_rank, geneticOut_);
		geneticOut_ << "All fitness: " << endl;
		GAPrinting::print_vector(fitness_energies, geneticOut_);		
	}
}


}