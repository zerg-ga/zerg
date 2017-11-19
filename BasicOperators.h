#ifndef BASICOPERATORS_H
#define BASICOPERATORS_H

#include <vector>
#include <fstream>

#include "Population.h"

namespace zerg{
class BasicOperators : public zerg::Population
{
public:
	BasicOperators(){};
	BasicOperators(int pop_size, int number_parameters);
	~BasicOperators();

	void startSimpleTest();

	//GA
	double getfitness(int ind_i);
	virtual void local_optimization(int ind_i);
	virtual bool create_individual(int creation_type, int target, int parent1, int parent2);
	virtual bool check_similarity(int target);
	inline int get_number_of_creation_methods(){return number_of_creation_methods;}
	virtual bool operatorAdministration(int method, const std::vector<double> &operatorPerformance);// modify operators
	virtual bool checkMinimum(int ind_i) { return true; }

protected:
	//data
	int number_of_creation_methods;
	int n_parameters;
	std::vector< std::vector<double> > x_vec;
	std::vector<double> energy;
	std::vector<double> random_individual_range_min;
	std::vector<double> random_individual_range_max;
	double tol_similarity; // can be modified by object

	void startBasicOperators();//final adjustments

	// operators 
	std::vector<double> generate_random_individual(); //always creation method = 0
	void make_mutation(int target, int parent);
	double mutationValue;
	void make_crossover_mean(int target,int parent1, int parent2);
	double crossoverWeight;
	void make_crossover_2_points(int target,int parent1, int parent2);
	void make_crossover_probability(int target,int parent1, int parent2);
	double crossoverProbability;

private:
	//procedures
	virtual void optimize(int ind_i); // base first
	void select_2_points(int &point1, int &point2);

	//similarity euclidean
	bool read_file_to_check_similarity(int ind_i);
	void append_to_similarity(int ind_i);
	std::ofstream write_similar_;
	std::ifstream read_similar_;
};

}

#endif


