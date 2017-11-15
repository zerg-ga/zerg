//#define useFuzzy

#include "FuzzyAdministration.h"

#include <iostream>

using namespace std;

namespace zerg{
FuzzyAdministration::~FuzzyAdministration()
{

#ifdef useFuzzy
	delete engine; // it deletes others
#endif
};

void FuzzyAdministration::setFuzzyRules(double maxEnergyVariarion_in, double maxCreationRateVariation_in)
{
	maxEnergyVariation = maxEnergyVariarion_in;

	maxCreationRateVariation = maxCreationRateVariation_in;

#ifdef useFuzzy
	using namespace fl;
	engine = new Engine;
	engine->setName("fuzzyGa");

	fitnessVariation = new InputVariable;
	fitnessVariation->setEnabled(true);
	fitnessVariation->setName("Variation");
	fitnessVariation->setRange(0.0000e0, 6.0000e0);

	fitnessVariation->addTerm(new Triangle("Magicarp", 0.000, 0.500, 0.750));
	fitnessVariation->addTerm(new Triangle("Terrible", 0.500, 1.250, 1.750));
	fitnessVariation->addTerm(new Triangle("Bad", 1.250, 1.750, 2.000));
	fitnessVariation->addTerm(new Triangle("Zero", 1.750, 2.000, 2.250));//2 is zero
	fitnessVariation->addTerm(new Triangle("Good", 2.000, 2.250, 2.750));
	fitnessVariation->addTerm(new Triangle("Excelent", 2.250, 2.750, 3.250));
	fitnessVariation->addTerm(new Triangle("Outstanding", 2.750, 3.250, 4.000));

	engine->addInputVariable(fitnessVariation);

	creationRateControl = new OutputVariable;
	creationRateControl->setEnabled(true);
	creationRateControl->setName("Rate");
	creationRateControl->setRange(0.000, 6.000);
	creationRateControl->fuzzyOutput()->setAccumulation(new Maximum);
	creationRateControl->setDefuzzifier(new Centroid(200));
	creationRateControl->setDefaultValue(12345);
	creationRateControl->setLockPreviousOutputValue(false);
	creationRateControl->setLockOutputValueInRange(false);

	creationRateControl->addTerm(new Triangle("toKill", 0.200, 0.400, 0.600));
	creationRateControl->addTerm(new Triangle("toGreatDecrease", 0.4000, 0.6000, 0.800));
	creationRateControl->addTerm(new Triangle("toDecrease", 0.600, 0.800, 1.000));
	creationRateControl->addTerm(new Triangle("toKeep", 0.800, 0.950, 1.050)); //1 is zero
	creationRateControl->addTerm(new Triangle("toIncrease", 1.000, 1.200, 1.400));
	creationRateControl->addTerm(new Triangle("toGreatIncrease", 1.200, 1.400, 1.600));
	creationRateControl->addTerm(new Triangle("toDominate", 1.400, 1.600, 1.800));

	engine->addOutputVariable(creationRateControl);

	ruleBlock = new RuleBlock;
	ruleBlock->setEnabled(true);
	ruleBlock->setName("");
	ruleBlock->setConjunction(fl::null);
	ruleBlock->setDisjunction(fl::null);
	ruleBlock->setActivation(new Minimum);

	ruleBlock->addRule(fl::Rule::parse("if Variation is Magicarp     then Rate is toKill", engine));
	ruleBlock->addRule(fl::Rule::parse("if Variation is Terrible     then Rate is toGreatDecrease", engine));
	ruleBlock->addRule(fl::Rule::parse("if Variation is Bad          then Rate is toDecrease", engine));
	ruleBlock->addRule(fl::Rule::parse("if Variation is Zero         then Rate is toKeep", engine));
	ruleBlock->addRule(fl::Rule::parse("if Variation is Good         then Rate is toIncrease", engine));
	ruleBlock->addRule(fl::Rule::parse("if Variation is Excelent     then Rate is toGreatIncrease", engine));
	ruleBlock->addRule(fl::Rule::parse("if Variation is Outstanding  then Rate is toDominate", engine));

	engine->addRuleBlock(ruleBlock);

	engine->configure("Minimum", "Maximum", "AlgebraicProduct", "AlgebraicSum", "Centroid");
#endif
}

double FuzzyAdministration::getCreateRateVariation(double methodMean)
{
	if(methodMean <= -maxEnergyVariation)
		return maxCreationRateVariation;
	else if(methodMean >= maxEnergyVariation)
		return -maxCreationRateVariation;
#ifndef useFuzzy
	double inclination = -maxCreationRateVariation / maxEnergyVariation;
	return inclination*(methodMean + maxEnergyVariation) + maxCreationRateVariation;

#else
	use namespace fl;
	scalar in = -methodMean + 2.0e0;
	fitnessVariation->setInputValue(in);
	engine->process();
	creationRateControl->defuzzify();
	return (creationRateControl->getOutputValue() - 1.0e0);
#endif
}

}