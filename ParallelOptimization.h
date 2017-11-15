#ifndef PARALLELOPTIMIZATION_H
#define PARALLELOPTIMIZATION_H

#include <vector>

#include "Population.h"

// fitness is evaluated here, the word parallel
// is just pressure
namespace zerg{
class ParallelOptimization
{
	int n_process;
public:
	void do_optimization(zerg::Population &pop, std::vector<int> &dead_individuals);
	inline void set_number_of_processor(int n_process_in){n_process=n_process_in;}

};

}

#endif

