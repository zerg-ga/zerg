#ifndef CREATION_H
#define CREATION_H

#include <vector>
#include <fstream>

#include "Population.h"
#include "Predator.h"
#include "ParallelOptimization.h"
#include "AdministrateCreation.h"
#include "StructOptions.h"
#include "Printing.h"

// porcentagem
// criacao pode ser assim:
// o ultimo perde 10 pontos e o penultimo perde 1.
// normaliza. Se o cara ficar de fora
// ele ganha um ponto.

namespace zerg{
class Creation
{
public:
	void initialize_creation(
		int pop_size,
		int number_creation_methods, 
		int n_process, 
		std::ofstream &geneticOut_,
		zerg::Printing * pPrinting_in,
		zerg::GaOptions &gaoptions,
		zerg::GaParameters &gaParam);

	void make_new_individuals(Population &pop, Predator &pred);

private:
	int number_methods;
	int insistOnSimilar;

	// creation methods administration
	std::vector<std::vector<int> > creation_methods; // quantity[0] and individuals to create
	std::vector<double> creation_rate;// how much will be created %

	void set_creation_methods(Predator &pred);
	void creation_from_methods(Population &pop, Predator &pred);
	void select_parents(Predator &pred, int &parent1, int &parent2);
	bool is_dead(Predator &pred,int parent);

	//printing functions
	void printSimilarityProblem(int method);

	zerg::ParallelOptimization go_parallel_;
	std::ofstream * pgeneticOut_;
	zerg::GaOptions * pgaoptions_;
	zerg::AdministrateCreation admin_;
	zerg::Printing * pPrinting_;

};

}

#endif

