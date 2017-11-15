#ifndef FUZZYGA_H
#define FUZZYGA_H

//#include <fl/Headers.h>

namespace zerg{
class FuzzyAdministration
{
public:
	FuzzyAdministration(){};
	~FuzzyAdministration();

	void setFuzzyRules(double maxEnergyVariation_in, double maxCrationRateVariation_in);

	double getCreateRateVariation(double methodMean);


private:
	double maxEnergyVariation;

	double maxCreationRateVariation;

	/*
	fl::Engine * engine;
	fl::InputVariable* fitnessVariation;
	fl::OutputVariable* creationRateControl;
	fl::RuleBlock* ruleBlock;
	*/
};

}

#endif



