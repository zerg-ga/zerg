#include "ClustersFitness.h"

#include "Fitness.h"
#include "WriteQuantumInput.h"
#include "../StructOptions.h"
#include "../Printing.h"
#include "../Random.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

using namespace std;
using namespace zerg;

ClustersFitness::ClustersFitness(
	Random * rand_in,
	GaParameters & gaParam,
	std::vector< std::string > &options_in,
	std::string gamessPath_in,
	std::string gamessScr_in,
	std::string nProc_in,
	Printing * pPrinting_in)
:ClustersOperators(
	gaParam.pop_size,
	gaParam.numberOfParameters,
	(int) gaParam.initialCreationRate.size(),
	pPrinting_in)
{
	options = options_in;
	optionsGamess = options_in;
	generationToChangeInteraction = gaParam.generationToChangeInteraction;
	iRestart = 0;
	interactionChanged = false;
	numberOfLocalMinimizations = 0;
	atomLabels = gaParam.atomLabels;
	setExperimentConditions(gaParam); //SBQT - EXPERIMENTS
	interactionType = gaParam.interactionPotentialType;
	if(interactionType == 4)
		changeInteraction(2);

	interactionParameters = gaParam.potentialParams;
	if (gaParam.restart)
	{
		readRestartFile();
		restartMax = restartEnergies.size();
	}
	else
	{
		remove("restart.ga");
		restartMax = 0;
	}

	startClustersOperators(rand_in, gaParam);

	gamessPath = gamessPath_in;
	gamessScr = gamessScr_in;
	nProc = nProc_in;
	scdo = gaParam.scdo;
	alfaMinGcdo = gaParam.alfaMinGcdo;
	alfaMaxGcdo = gaParam.alfaMaxGcdo;
	wGcdo = gaParam.wGcdo;
	tetaMinTwisto = gaParam.tetaMinTwisto;
	tetaMaxTwisto = gaParam.tetaMaxTwisto;
	contractionMinMtco = gaParam.contractionMinMtco;
	contractionMaxMtco = gaParam.contractionMaxMtco;

	seed = gaParam.seed;
	experimentMethod = gaParam.experimentMethod;
	bool aux;
	if(sim_.printLevel > 1)
		sim_.printLevel = 1;
	sim_.printInitialTitle();
	pPrinting_->startingFirstIndividual();
	aux = create_individual(0, 0, 0, 0); //method 0 always random

	local_optimization(0);

	pPrinting_->endOfFirstIndividual();

	bool addInd;
	for(int i=1; i<gaParam.pop_size; i++)
	{
		addInd = false;
		for (int k = 0; k < gaParam.insistOnSimilar; k++)
		{
			aux = create_individual(0, i, 0, 0); //method 0 always random
			if (!checkInitialSimilarity(i))
				break;
		}
		local_optimization(i);
		if (checkInitialSimilarity(i))
			energy[i] = 1.0e99;
	}
	pPrinting_->endOfInitialPopulation();
	sim_.printLevel = gaParam.similarityDebugLevel;
	sim_.printTitle();
}

ClustersFitness::~ClustersFitness(){}


void ClustersFitness::local_optimization(int ind_i)
{
	if (iRestart < restartMax)
	{
		energy[ind_i] = restartEnergies[iRestart];
		x_vec[ind_i] = restartCoordinates[iRestart];
		iRestart++;
	}
	else
	{
		if((generation == generationToChangeInteraction) && (!interactionChanged))
		{
			changeInteraction(3);
			interactionChanged = true;
			for(size_t i = 0; i < energy.size(); i++)
			{
				optimize(i);
				saveIndividual(i);
			}
		}

		optimize(ind_i);
		saveIndividual(ind_i);
	}
	translateToGeometricCenter(x_vec[ind_i]);
	appendTosimilarity(ind_i);
}

void ClustersFitness::optimize(int ind_i)
{
	// you have:
	// x_vec[ind_i] is a vector<double> -> to be optimized
	// I want:
	// energy[ind_i] -> fitness function

	sim_.printNewBfgsInd();

	Fitness fit_;
	if (options.size() == 0)
	{
		energy[ind_i] = fit_.optimizeEmpiricalPotential(
			x_vec[ind_i],
			interactionType,
			interactionParameters,
			atomTypes,
			&sim_);
	}
	else
	{
		if(options[0] == "gamess")
		{
			energy[ind_i] = fit_.runGamess(
				x_vec[ind_i],
				options,
				gamessPath,
				gamessScr,
				nProc,
				&sim_);
		}
		else if(options[0] == "mopac")
		{
			energy[ind_i] = fit_.runMopac(
                                x_vec[ind_i],
                                options,
				atomLabels);
		}
		else
		{
			cout << "Quantum method not found - exiting" << endl;
			exit(1);
		}


	}

	sim_.printEndIntoBfgs();

	numberOfLocalMinimizations++;

	if (makeExperiment)
		endExperimentConditions(ind_i);
}

void ClustersFitness::printAllIndividuals(string fileName)
{
	// ordering
	int pop_size = energy.size();
	vector<double> fitness_energies(pop_size);
	vector<int> fitness_rank(pop_size);
	for (int i = 0; i<pop_size; i++)
	{
		fitness_energies[i] = energy[i];
		fitness_rank[i] = i;
	}
	vector<int> vector_order = auxMath_.vector_ordering(fitness_energies);
	auxMath_.vector_ordering_with_instructions(fitness_rank, vector_order);

	ofstream printAll_(fileName.c_str());
	int nAtoms = x_vec[0].size() / 3;
	for (size_t ind = 0; ind < x_vec.size(); ind++)
	{
		int best = fitness_rank[ind];
		printAll_ << nAtoms << endl << setprecision(16) << energy[best] << endl;
		for (int i = 0; i < nAtoms; i++)
		{
			printAll_ << atomLabels[i] <<  "  ";

			printAll_ << x_vec[best][i] << "  "
				<< x_vec[best][i + nAtoms] << "  "
				<< x_vec[best][i + 2 * nAtoms] << "  "
				<< endl;
		}
	}
	printAll_.close();
}

void ClustersFitness::saveIndividual(int ind_i)
{
	saveIndividual_.open("restart.ga", std::ofstream::out | std::ofstream::app);
	int nAtoms = x_vec[0].size() / 3;
	saveIndividual_
		<< nAtoms << endl
		<< setprecision(16) << energy[ind_i] << endl;
	for (int i = 0; i < nAtoms; i++)
	{
		saveIndividual_ << "H "
		<< x_vec[ind_i][i] << "  "
		<< x_vec[ind_i][i + nAtoms] << "  "
		<< x_vec[ind_i][i + 2 * nAtoms] << "  "
		<< endl;
	}
	saveIndividual_.close();
}

void ClustersFitness::readRestartFile()
{
	ifstream restart_("restart.ga");
	string auxline;
	string dummy;
	int nAtoms;
	while (getline(restart_, auxline))
	{
		stringstream line0;
		line0 << auxline;
		line0 >> nAtoms;
		vector<double> coord(3 * nAtoms);
		double energyI;
		getline(restart_, auxline);
		stringstream lineEnerg;
		lineEnerg << auxline;
		lineEnerg >> energyI;
		for (int i = 0; i < nAtoms; i++)
		{
			getline(restart_, auxline);
			stringstream line;
			line << auxline;
			line >> dummy
				>> coord[i]
				>> coord[i + nAtoms]
				>> coord[i + 2 * nAtoms];
		}
		if (coord.size() != nAtoms)
		{
			restartCoordinates.push_back(coord);
			restartEnergies.push_back(energyI);
		}
	}
	restart_.close();
}

void ClustersFitness::setExperimentConditions(GaParameters &gaParam)
{
	makeExperiment = gaParam.isToMakeExperiment;
	globalMinimaEnergy = gaParam.experimentFinalEnergy;
	maxMinimizations = gaParam.experimentMaxSteps;
}

void ClustersFitness::endExperimentConditions(int ind_i)
{
	double energyI = energy[ind_i];
	if (energyI <= globalMinimaEnergy)
	{
		ofstream fileCsv_;
		fileCsv_.open(("cluster-result-" + experimentMethod + ".csv").c_str(), std::ofstream::out | std::ofstream::app);
		fileCsv_ << experimentMethod << " ; " << seed << " ; " << numberOfLocalMinimizations << endl;
		fileCsv_.close();
		saveIndividual(ind_i);
		exit(0);
	}
	else if (numberOfLocalMinimizations > maxMinimizations)
	{
		ofstream fileCsv_;
		fileCsv_.open(("cluster-result-" + experimentMethod + ".csv").c_str(), std::ofstream::out | std::ofstream::app);
		fileCsv_ << experimentMethod << " ; " << seed << " ; " << numberOfLocalMinimizations << endl;
		fileCsv_.close();
		exit(0);
	}
}

double ClustersFitness::checkMinimum(int ind_i)
{
	if (options.size() != 0)
	{
		if(options[0] == "gamess")
		{
			Fitness fit_;
			energy[ind_i] = fit_.runGamess(
				x_vec[ind_i],
				options,
				gamessPath,
				gamessScr,
				nProc);
	
			double frequency = fit_.runGamessFrequency(
				numberOfLocalMinimizations,
				x_vec[ind_i],
				options,
				gamessPath,
				gamessScr,
				nProc);

			if (frequency < 0.0e0)
				energy[ind_i] = 1.0e99;

			return frequency;
		}
		else
		{
			return 1.0e0;

		}
	}
	return 1.0e0;
}


void ClustersFitness::changeInteraction(int newType)
{
	switch(newType)
	{
		case 0:
		case 1:
			options.clear();
			interactionType = newType;
			break;
		case 2:
			options.clear();
	                options.push_back("mopac");
        	        options.push_back(optionsGamess[1]);
			options.push_back("PM7");
			options.push_back("");
			options.push_back("");
			break;
		case 3:
			options = optionsGamess;
			break;
		default:
			cout << "new interaction type not found - exiting" << endl;
			exit(1);
	}
}

void ClustersFitness::printBfgsSteps()
{
	sim_.printBfgsSteps();
}
