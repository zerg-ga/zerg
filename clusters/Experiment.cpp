#include "Experiment.h"

#include <iostream>

#include "ReadGaInput.h"
#include "ClustersFitness.h"
#include "../AuxMathGa.h"
#include "../AuxMath.h"
#include "../GeneticAlgorithm.h"
#include "../Printing.h"

using namespace zerg;
using namespace std;

/*
objeto que seta os parametros e roda o ga umas
5 vezes com os mesmos parametros mas com sementes diferentes.

*/

Experiment::Experiment() {}

Experiment::~Experiment() {}

void Experiment::makeExperiment(int seed, string experimentMethod, vector<double> & additionalParams)
{
	Printing printing_;

	int nAtoms = 26;

	ReadGaInput readGa_ (&printing_);
		
	readGa_.setExperimentDefaults(seed);

	zerg::GaParameters gaParam = readGa_.getGaParameters();

	// best parameters
	AuxMath auxMath_;
	gaParam.tetaMinTwisto = auxMath_._pi * 0.1e0;
	gaParam.tetaMaxTwisto = auxMath_._pi * 0.5e0;
	gaParam.alfaMinGcdo = 0.2;
	gaParam.alfaMaxGcdo = 0.7;
	gaParam.wGcdo = 2;

	if (experimentMethod == "AutoAdjust")
	{
		gaParam.adminMaxCreationVariation = additionalParams[0];
		gaParam.adminLargeEnergyVariation = additionalParams[1];
	}
	else if (experimentMethod == "TwistOperator")
	{
		AuxMath auxMath_;
		gaParam.tetaMinTwisto = auxMath_._pi * additionalParams[0];
		gaParam.tetaMaxTwisto = auxMath_._pi * additionalParams[1];
	}
	else if (experimentMethod == "GeometricCenterDisplacement")
	{
		gaParam.alfaMinGcdo = additionalParams[0];
		gaParam.alfaMaxGcdo = additionalParams[1];
		gaParam.wGcdo = additionalParams[2];
	}

	gaParam.experimentMethod = experimentMethod;

	setExperiment(experimentMethod, gaParam);

	AuxMathGa::set_seed(gaParam.seed);

	gaParam.numberOfParameters = 3 * nAtoms;

	vector<string> options = readGa_.getOptions();

	ClustersFitness clFit_(
		gaParam,
		options,
		readGa_.getGamessPath(),
		readGa_.getGamessScr(),
		readGa_.getGamessNprocess());

	clFit_.setExperimentConditions(-108.315e0, 3000);

	GeneticAlgorithm ga1(clFit_, gaParam, &printing_);

	ga1.ga_start();



	// avaliar a melhora da populacao como um todo.
	// nao apenas o melhor individuo
	// o 26 e um caso classico, acho que eu poderia me concentrar nele.
	// a qualidade dos operadores ao longo do tempo
	// se as variacoes sao melhores no inicio da simulacao, no final, ou nao depende e tal.
	// adicionar um Nlm e um Ncalls
}


void Experiment::setExperiment(std::string experimentMethod, zerg::GaParameters & gaParam)
{
	for (size_t i = 1; i < gaParam.initialCreationRate.size(); i++)
		gaParam.initialCreationRate[i] = 0.0e0;

	if (experimentMethod == "AutoAdjust")
	{
		for (size_t i = 1; i < gaParam.initialCreationRate.size(); i++)
			gaParam.initialCreationRate[i] = 1.0e0 / ((int)gaParam.initialCreationRate.size()-1);
	}
	else if (experimentMethod == "TwistOperator")
		gaParam.initialCreationRate[3] = 1.0e0;
	else if (experimentMethod == "DeavenHoCutSplice")
		gaParam.initialCreationRate[4] = 1.0e0;
	else if (experimentMethod == "SphereCrossover")
		gaParam.initialCreationRate[1] = 1.0e0;
	else if (experimentMethod == "GeometricCenterDisplacement")
		gaParam.initialCreationRate[2] = 1.0e0;
	else if (experimentMethod == "AngularOperator")
		gaParam.initialCreationRate[5] = 1.0e0;
	else if (experimentMethod == "AngularSurface")
		gaParam.initialCreationRate[6] = 1.0e0;
	else
	{
		cout << "A NUMERACAO DOS OPERADORES MUDOU - CUIDADO" << endl;
		exit(1);
	}
	return;





	if (experimentMethod == "Imigration")
		gaParam.initialCreationRate[0] = 0.5e0;
	else if (experimentMethod == "CrossoverMean")
		gaParam.initialCreationRate[1] = 0.5e0;
	else if (experimentMethod == "CrossoverTwoPoints")
		gaParam.initialCreationRate[2] = 0.5e0;
	else if (experimentMethod == "SphereCrossover")
		gaParam.initialCreationRate[5] = 0.5e0;
	else if (experimentMethod == "CartesianDisplacement")
		gaParam.initialCreationRate[6] = 0.5e0;
	else if (experimentMethod == "GeometricCenterDisplacement")
		gaParam.initialCreationRate[7] = 0.5e0;
	else if (experimentMethod == "TwistOperator")
		gaParam.initialCreationRate[8] = 0.5e0;
	else if (experimentMethod == "AngularOperator")
		gaParam.initialCreationRate[9] = 0.5e0;
	else if (experimentMethod == "AngularSurface")
		gaParam.initialCreationRate[10] = 0.5e0;
	else if (experimentMethod == "MoveToCenter")
		gaParam.initialCreationRate[11] = 0.5e0;
	else if (experimentMethod == "ModifiedAngularSurface")
		gaParam.initialCreationRate[12] = 0.5e0;
	else if (experimentMethod == "DeavenHoCutSplice")
		gaParam.initialCreationRate[13] = 0.5e0;
	else if (experimentMethod == "SmallAutoAdjust")
	{
		for (size_t i = 0; i < gaParam.initialCreationRate.size(); i++)
			gaParam.initialCreationRate[i] = 1.0e0 / (int)gaParam.initialCreationRate.size();
		gaParam.adminMaxCreationVariation = 0.2e0;
	}
	else if (experimentMethod == "MediumAutoAdjust")
	{
		for (size_t i = 0; i < gaParam.initialCreationRate.size(); i++)
			gaParam.initialCreationRate[i] = 1.0e0 / (int)gaParam.initialCreationRate.size();
		gaParam.adminMaxCreationVariation = 0.5e0;
	}
	else if (experimentMethod == "LargeAutoAdjust")
	{
		for (size_t i = 0; i < gaParam.initialCreationRate.size(); i++)
			gaParam.initialCreationRate[i] = 1.0e0 / (int)gaParam.initialCreationRate.size();
		gaParam.adminMaxCreationVariation = 0.9e0;
	}
	else if (experimentMethod == "ImigrationOnly")
	{
		gaParam.initialCreationRate[0] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "CrossoverMeanOnly")
	{
		gaParam.initialCreationRate[1] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "CrossoverTwoPointsOnly")
	{
		gaParam.initialCreationRate[2] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "SphereCrossoverOnly")
	{
		gaParam.initialCreationRate[5] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "CartesianDisplacementOnly")
	{
		gaParam.initialCreationRate[6] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "GeometricCenterDisplacementOnly")
	{
		gaParam.initialCreationRate[7] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "TwistOperatorOnly")
	{
		gaParam.initialCreationRate[8] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "AngularOperatorOnly")
	{
		gaParam.initialCreationRate[9] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "AngularSurfaceOnly")
	{
		gaParam.initialCreationRate[10] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "MoveToCenterOnly")
	{
		gaParam.initialCreationRate[11] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "ModifiedAngularSurfaceOnly")
	{
		gaParam.initialCreationRate[12] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "DeavenHoCutSpliceOnly")
	{
		gaParam.initialCreationRate[13] = 1.0e0;
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "Mutation")
	{
		gaParam.initialCreationRate[3] = 0.7e0;
	}
	else if (experimentMethod == "CrossoverProbability")
	{
		gaParam.initialCreationRate[4] = 0.8e0;
	}
	else if (experimentMethod == "MutationOnly")
	{
		gaParam.initialCreationRate[3] = 1.0e0;
		gaParam.initialCreationRate[4] = 0.0e0;
	}
	else if (experimentMethod == "CrossoverProbabilityOnly")
	{
		gaParam.initialCreationRate[3] = 0.0e0;
		gaParam.initialCreationRate[4] = 1.0e0;
	}
	else
	{
		cout << "method not found" << endl;
		exit(1);
	}
}
