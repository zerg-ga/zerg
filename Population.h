#ifndef POPULATION_H
#define POPULATION_H

#include <vector>

namespace zerg{
class Population
{
public:
	virtual ~Population(){};

	virtual double getfitness(int ind_i) = 0;

	virtual void local_optimization(int ind_i) = 0;// fitness

	virtual bool create_individual(int creation_type, int target, int parent1, int parent2) = 0; // all creations

	virtual int get_number_of_creation_methods() = 0;

	virtual bool check_similarity(int target) = 0;

	virtual bool operatorAdministration(int method, const std::vector<double> &operatorPerformance) = 0;// modify operators

};

}
#endif


