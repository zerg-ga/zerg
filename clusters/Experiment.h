#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <string>
#include "ReadGaInput.h"

#include "../StructOptions.h"

class Experiment
{
public:
	Experiment();

	~Experiment();

	void makeExperiment(int seed, std::string experimentMethod, std::vector<double> & additionalParams);

private:

	void setExperiment(std::string experimentMethod, zerg::GaParameters & gaParam);

};

#endif

