#include "ReadGaInput.h"


#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdlib.h>

#include "../AuxMath.h"
#include "../Printing.h"
#include "WriteQuantumInput.h"

using namespace zerg;
using namespace std;

ReadGaInput::~ReadGaInput(){}

ReadGaInput::ReadGaInput(Printing * pPrinting_in)
{
	inputName = "GaInput.txt";
	pPrinting_ = pPrinting_in;
}

void ReadGaInput::readGaInput()
{
	setDefaults();
	string atom1Default = "H";
	string atom2Default = "N";
	string atom3Default = "C";

	ifstream input_(inputName.c_str());
	vector<string> baseFiles;
	string auxline, type, equal, value, projectName, gamessHeader;
	string interactionPotential;

	pPrinting_->gaStartInputReading();

	while (getline(input_, auxline))
	{
		pPrinting_->showInputLines(auxline);
		stringstream line;
		line << auxline;
		line >> type >> equal >> value;
		if (type == "%")
			continue;
		if (equal != "=")
		{
			cout << "wrong input format" << endl;
			exit(1);
		}
		stringstream convert;
		convert << value;
		if (type == "seed")
			convert >> gaParam.seed;
		else if (type == "restart")
			gaParam.restart = convertInputToBool(convert);
		else if (type == "number_of_cores")
			convert >> gamessNproc;
		else if (type == "project_name")
			convert >> projectName;
		else if (type == "printing_debug")
			convert >> gaParam.printingDebug;
		else if (type == "population_size")
		{
			convert >> gaParam.pop_size;
			if (gaParam.pop_size < 4)
			{
				cout << "Sorry, population_size must be higher than four" << endl 
					<< "Exiting" << endl;
				exit(1);
			}
		}
		else if (type == "maximum_number_of_generations")
			convert >> gaParam.maxGeneration;
		else if (type == "highlander_initial_energy")
			convert >> gaParam.highlanderInitialFitness;
		else if (type == "highlander_max_iteration")
			convert >> gaParam.highlanderMaxIteration;
		else if (type == "highlander_threshold")
			convert >> gaParam.highlanderThreshold;
		else if (type == "large_energy_value")
			convert >> gaParam.adminLargeEnergyVariation;
		else if (type == "maximum_creation_operators_boost")
			convert >> gaParam.adminMaxCreationVariation;
		else if (type == "predator_method")
			convert >> gaParam.predatorMethod;
		else if (type == "number_of_predated")
			convert >> gaParam.numberOfKilling;
		else if (type == "mutation_value")
			convert >> gaParam.mutationValue;
		else if (type == "crossover_weight")
			convert >> gaParam.crossoverWeight;
		else if (type == "crossover_probability")
			convert >> gaParam.corssoverProbability;
		else if (type == "sCartesianDisplacementOperator")
			convert >> gaParam.scdo;
		else if (type == "alfaMinGeometricDisplacement")
			convert >> gaParam.alfaMinGcdo;
		else if (type == "alfaMaxGeometricDisplacement")
			convert >> gaParam.alfaMaxGcdo;
		else if (type == "wGeometricDisplacement")
			convert >> gaParam.wGcdo;
		else if (type == "tetaMinTwistOperator")
			convert >> gaParam.tetaMinTwisto;
		else if (type == "tetaMaxTwistOperator")
			convert >> gaParam.tetaMaxTwisto;
		else if (type == "contractionMinMoveToCenter")
			convert >> gaParam.contractionMinMtco;
		else if (type == "checkMinimum")
			gaParam.isToCheckMinimum = convertInputToBool(convert);
		else if (type == "contractionMaxMoveToCenter")
			convert >> gaParam.contractionMaxMtco;
		else if (type == "activateIntoBfgs")
			gaParam.activateIntoBfgs = convertInputToBool(convert);
		else if (type == "similarityMethod")
			convert >> gaParam.similarityMethod;
		else if (type == "removeSimilarStructures")
			convert >> gaParam.removeSimilarStructures;
		else if (type == "similarityDebugLevel")
			convert >> gaParam.similarityDebugLevel;
		else if (type == "tolSimilarity")
			convert >> gaParam.tolSimilarity;
		else if (type == "bestIndividualSize")
			convert >> gaParam.bestIndividualSize;
		else if (type == "energyReturnBfgs")
			convert >> gaParam.energyReturnBfgs;
		else if (type == "number_of_atoms")
		{
			convert >> gaParam.numberOfParameters;
			gaParam.numberOfParameters *= 3;
			gaParam.nAtomTypes1 = gaParam.numberOfParameters / 3;
			gaParam.nAtomTypes2 = 0;
			gaParam.nAtomTypes3 = 0;
		}
		else if (type == "n_atom_type_1")
			convert >> gaParam.nAtomTypes1;
		else if (type == "n_label_type_1")
			convert >> atom1Default;
		else if (type == "n_atom_type_2")
			convert >> gaParam.nAtomTypes2;
		else if (type == "n_label_type_2")
			convert >> atom2Default;
		else if (type == "n_atom_type_3")
			convert >> gaParam.nAtomTypes3;
		else if (type == "n_label_type_3")
			convert >> atom3Default;
		else if(type == "user_defined_method")
		{
			int userMethod;
			convert >> userMethod;
			userDefinedSet(userMethod);
		}
		else if (type == "isToMakeExperiment")
			gaParam.isToMakeExperiment = convertInputToBool(convert);
		else if (type == "experimentFinalEnergy")
			convert >> gaParam.experimentFinalEnergy;
		else if (type == "experimentMaxSteps")
			convert >> gaParam.experimentMaxSteps;
		else if (type == "gamma_creation_radius")
			convert >> gaParam.gammaInitializeAtoms;
		else if (type == "gamma_creation_radius")
			convert >> gaParam.gammaInitializeAtoms;
		else if (type == "gamma_creation_radius")
			convert >> gaParam.gammaInitializeAtoms;
		else if (type == "radius_factor")
			convert >> gaParam.rcaInitializeAtoms;
		else if (type == "max_distance_between_atoms")
			convert >> gaParam.maxDistance;
		else if (type == "min_distance_between_atoms")
			convert >> gaParam.minDistance;
		else if (type == "iterations_to_repeat_if_it_is_similar")
			convert >> gaParam.insistOnSimilar;
		else if (type == "interaction_potential")
			convert >> interactionPotential;
		else if (type == "generation_change_interaction")
			convert >> gaParam.generationToChangeInteraction;
		else if (type == "number_of_parameters")
		{
			int nParameters;
			convert >> nParameters;
			for (int i = 0; i < nParameters; i++)
			{
				string parametersLine, typeParam;
				double paramValue;
				getline(input_, parametersLine);
				stringstream lineBase;
				lineBase << parametersLine;
				lineBase >> typeParam >> equal >> paramValue;

				if (typeParam == "%")
					continue;
				else if (equal != "=")
				{
					cout << "wrong input format" << endl;
					exit(1);
				}
				gaParam.potentialParams.push_back(paramValue);
			}
		}
		else if (type == "gamess_executable_path")
			convert >> gamessPath;
		else if (type == "gamess_scr_path")
			convert >> gamessScr;
		else if (type == "bases_files_number")
		{
			int nBases;
			convert >> nBases;
			for (int i = 0; i < nBases; i++)
			{
				string auxBasesName;
				string typeBase, baseName;
				getline(input_, auxBasesName);
				stringstream lineBase;
				lineBase << auxBasesName;
				lineBase >> typeBase >> equal >> baseName;
				if (typeBase == "%")
					continue;
				else if (equal != "=")
				{
					cout << "wrong input format" << endl;
					exit(1);
				}
				else if (typeBase == "base")
					baseFiles.push_back(baseName);
			}
		}
		else if (type == "gamess_header_file")
			convert >> gamessHeader;
		else
		{
			cout << type << "  command not found, check manual for details " << endl;
			exit(1);
		}
	}
	input_.close();

	if(gaParam.numberOfKilling == 0)
		gaParam.numberOfKilling = (int)gaParam.pop_size / 4;
	if (gaParam.numberOfKilling >= gaParam.pop_size)
	{
		cout << "number_of_predated must be lower than population_size" << endl
			<< "exiting" << endl;
		exit(1);
	}
	if (gaParam.numberOfKilling < 1)
	{
		cout << "number_of_predated must be equal to one or higher" << endl
			<< "exiting" << endl;
		exit(1);
	}

	//show final options
	pPrinting_->showAllParameters(
		gaParam,
		gamessNproc,
		projectName,
		interactionPotential,
		gamessPath,
		gamessScr,
		gamessHeader,
		baseFiles);

	int nAtoms = gaParam.numberOfParameters / 3;
	int nTotalTypes = gaParam.nAtomTypes1 + gaParam.nAtomTypes2 + gaParam.nAtomTypes3;
	if (nAtoms != nTotalTypes)
	{
		cout << "total entered atoms: " << nAtoms << endl
			<< "sum of atoms on atom types: " << nTotalTypes << endl
			<< "values don't match, check input for errors" << endl;
		exit(1);
	}

	vector<string> newLabelVec;
	gaParam.atomLabels = newLabelVec;
	vector<int> atomTypesTemp(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		if (i < gaParam.nAtomTypes1)
		{
			atomTypesTemp[i] = 0;
			gaParam.atomLabels.push_back(atom1Default);
		}
		else if (i < (gaParam.nAtomTypes1 + gaParam.nAtomTypes2))
		{
			atomTypesTemp[i] = 1;
			gaParam.atomLabels.push_back(atom2Default);
		}
		else
		{
			atomTypesTemp[i] = 2;
			gaParam.atomLabels.push_back(atom3Default);
		}
	}
	gaParam.atomTypes = atomTypesTemp;

	if (interactionPotential == "lennardjones")
		gaParam.interactionPotentialType = 0;
	else if (interactionPotential == "gupta")
		gaParam.interactionPotentialType = 1;
	else if (interactionPotential == "mopac")
	{
		gaParam.interactionPotentialType = 2;
		options.push_back("mopac");
		options.push_back(projectName);
		options.push_back("PM7");
		options.push_back("");
		options.push_back("");
	}
	else if ((interactionPotential == "gamess") || (interactionPotential == "MopacGamess"))
	{
		gaParam.interactionPotentialType = 3;
		if(interactionPotential == "MopacGamess")
			gaParam.interactionPotentialType = 4;

		options.push_back("gamess");
		options.push_back(projectName);
		ifstream gamHeader_(("auxFiles/" + gamessHeader).c_str());
		while (getline(gamHeader_, auxline))
		{
			if (auxline.find("EndOfGamessHeader") != string::npos)
				break;
			options.push_back(auxline);
		}
		gamHeader_.close();
		options.push_back("EndOfHeader");
		bool haveOptimizeFlag = false;
		for (size_t i = 0; i < options.size(); i++)
		{
			if (options[i].find("RUNTYP=OPTIMIZE") != string::npos)
				haveOptimizeFlag = true;
		}
		if (!haveOptimizeFlag)
		{
			cout << "flag RUNTYP=OPTMIZE dont found." << endl
				<< "it must be uppercase, check header and input for errors" << endl;
			exit(1);
		}
	
		baseFilesOrdering(
			baseFiles, 
			gaParam.nAtomTypes1, 
			gaParam.nAtomTypes2, 
			gaParam.nAtomTypes3);

		for (size_t k = 0; k < baseFiles.size(); k++)
			options.push_back(baseFiles[k]);
		options.push_back("EndOfBasis");
		options.push_back("NoECP");

		options.push_back("End of options");
		options.push_back("Header name");
		options.push_back(gamessHeader);

		vector<int> auxAtomTypes(nAtoms);
		gaParam.atomLabels.clear();
		WriteQuantumInput writeInp_(options);
		for (int i = 0; i < nAtoms; i++)
		{
			if (i < gaParam.nAtomTypes1)
			{
				auxAtomTypes[i] = 0;
			}
			else if (i < (gaParam.nAtomTypes1 + gaParam.nAtomTypes2))
			{
				auxAtomTypes[i] = 1;
			}
			else
			{
				auxAtomTypes[i] = 2;
			}
			gaParam.atomLabels.push_back(writeInp_.getAtomName(i));
		}
		gaParam.atomTypes = auxAtomTypes;

		pPrinting_->endOfGamessOptions();
	}
	else
	{
		cout << "Interaction potential not found - exiting" << endl;
		exit(1);
	}
}

bool ReadGaInput::convertInputToBool(std::stringstream &convert)
{
	int readActivate;
	convert >> readActivate;
	return readActivate != 0;
}

void ReadGaInput::baseFilesOrdering(
	std::vector<std::string> & baseFiles,
	int & n1,
	int & n2,
	int & n3)
{

	string firstBaseFile = baseFiles[0];
	bool tradePositions = true;	
	while(tradePositions)
	{
		tradePositions = false;	
		for(size_t i = 0; i < baseFiles.size() - 1; i++)
		{
			if((baseFiles[i + 1] == firstBaseFile) &&
				(baseFiles[i] != firstBaseFile))
			{
				tradePositions = true;
				string auxBase = baseFiles[i];
				baseFiles[i] = baseFiles[i + 1];
				baseFiles[i + 1] = auxBase;
			}
		}
	}

	string lastBaseFile = baseFiles[baseFiles.size() - 1];
	if(lastBaseFile == firstBaseFile)
		tradePositions = false;	
	else
		tradePositions = true;	

	while(tradePositions)
	{
		tradePositions = false;	
		for(size_t i = 0; i < baseFiles.size() - 1; i++)
		{
			if(baseFiles[i] == firstBaseFile)
				continue;
			else if(baseFiles[i] == lastBaseFile)
				continue;
			else
			{
				tradePositions = true;
				string auxBase = baseFiles[i];
				baseFiles[i] = baseFiles[i + 1];
				baseFiles[i + 1] = auxBase;
			}
		}
	}

	n1 = 0;
	n2 = 0;
	n3 = 0;
	for(size_t i = 0; i < baseFiles.size(); i++)
	{
		if(baseFiles[i] == firstBaseFile)
			n1++;	
		else if(baseFiles[i] == lastBaseFile)
			n2++;
		else
			n3++;	
	}
}

void ReadGaInput::setInputInformations(GaParameters gaParam_in)
{
	gaParam = gaParam_in;
}

vector<string> ReadGaInput::getOptions()
{
	return options;
}

void ReadGaInput::userDefinedSet(int userMethod)
{
	if(userMethod == 17) {  // auto5
		int nOperators = 5;
		double nOperatorsRate = 1.0e0 / (double)nOperators;
		gaParam.initialCreationRate.resize(nOperators);
		for(int i = 0; i < nOperators; i++)
			gaParam.initialCreationRate[i] = nOperatorsRate;
		return;
	} else if(userMethod == 18) {  // auto7
		int nOperators = 7;
		double nOperatorsRate = 1.0e0 / (double)nOperators;
		gaParam.initialCreationRate.resize(nOperators);
		for(int i = 0; i < nOperators; i++)
			gaParam.initialCreationRate[i] = nOperatorsRate;
		return;
	} else  if(userMethod == 19) {  // auto13
		int nOperators = 13;
		double nOperatorsRate = 1.0e0 / (double)nOperators;
		gaParam.initialCreationRate.resize(nOperators);
		for(int i = 0; i < nOperators; i++)
			gaParam.initialCreationRate[i] = nOperatorsRate;
		return;
	}

	int EXCHANGE_OPERATOR = 13;
	gaParam.initialCreationRate.resize(EXCHANGE_OPERATOR + 1);
	gaParam.adminMaxCreationVariation = 0.0e0;
	gaParam.adminLargeEnergyVariation = 0.0e0;
	for(size_t i = 0; i < gaParam.initialCreationRate.size(); i++)
		gaParam.initialCreationRate[i] = 0.0e0;
	switch(userMethod)
	{
		case 1:
			// Deaven and Ho genetic algorithm
        		gaParam.initialCreationRate[0] = 0.1e0;
        		gaParam.initialCreationRate[EXCHANGE_OPERATOR] = 0.1e0;
        		gaParam.initialCreationRate[4] = 0.7e0;
        		gaParam.initialCreationRate[7] = 0.1e0;
			break;
		case 2:
        		gaParam.initialCreationRate[0] = 1.0e0;
			break;
		case 3:
        		gaParam.initialCreationRate[1] = 1.0e0;
			break;
		case 4:
        		gaParam.initialCreationRate[2] = 1.0e0;
			break;
		case 5:
        		gaParam.initialCreationRate[3] = 1.0e0;
			break;
		case 6:
        		gaParam.initialCreationRate[4] = 1.0e0;
			break;
		case 7:
        		gaParam.initialCreationRate[5] = 1.0e0;
			break;
		case 8:
        		gaParam.initialCreationRate[6] = 1.0e0;
			break;
		case 9:
        		gaParam.initialCreationRate[7] = 1.0e0;
			break;
		case 10:
        		gaParam.initialCreationRate[8] = 1.0e0;
			break;
		case 11:
        		gaParam.initialCreationRate[9] = 1.0e0;
			break;
		case 12:
        		gaParam.initialCreationRate[10] = 1.0e0;
			break;
		case 13:
        		gaParam.initialCreationRate[11] = 1.0e0;
			break;
		case 14:
        		gaParam.initialCreationRate[12] = 1.0e0;
			break;
		case 15:
        		gaParam.initialCreationRate[13] = 1.0e0;
			break;
		case 16: // Previous genetic algorithm build
        		gaParam.initialCreationRate[0] = 0.1e0;
        		gaParam.initialCreationRate[9] = 0.7e0;
        		gaParam.initialCreationRate[12] = 0.2e0;
			break;
		default:
			cout << "user method not found" << endl;
			exit(1);
	}
}

void ReadGaInput::setDefaults()
{
	gaParam.seed = 3;
	gaParam.restart = false;
	gaParam.pop_size = 40;
	gaParam.maxGeneration = 300;
	gaParam.highlanderInitialFitness = 1.0e99;
	gaParam.highlanderMaxIteration = 50;
	gaParam.highlanderThreshold = 1.0e-6;
	gaParam.insistOnSimilar = 30;
	gaParam.adminLargeEnergyVariation = 2.0e0;
	gaParam.adminMaxCreationVariation = 0.9e0;
	gaParam.predatorMethod = 0;
	gaParam.n_process = 1;
	gaParam.mutationValue = 0.1e0;
	gaParam.crossoverWeight = 0.7e0;
	gaParam.corssoverProbability = 0.7e0;
	gaParam.numberOfParameters = 15;
	gaParam.gammaInitializeAtoms = 0.2;
	gaParam.rcaInitializeAtoms = 1.0;
	gaParam.maxDistance = 1.0e99;
	gaParam.minDistance = 0.2e0;
	gaParam.printingDebug = 0;
	gaParam.isToCheckMinimum = false;

	AuxMath auxMath_;

	// OPERATORS PARAMETERS
	gaParam.alfaMinGcdo = 0.2e0;
	gaParam.alfaMaxGcdo = 0.7e0;
	gaParam.wGcdo = 2;
	gaParam.tetaMinTwisto = 0.1e0 * auxMath_._pi;
	gaParam.tetaMaxTwisto = 0.5e0 * auxMath_._pi;
	gaParam.scdo = 0.2e0;
	gaParam.contractionMinMtco = 0.1e0;
	gaParam.contractionMaxMtco = 0.8e0;
	gaParam.numberOfKilling = 0;

	// SIMILARITY PARAMETERS
	gaParam.activateIntoBfgs = false;
	gaParam.similarityMethod = -1;
	gaParam.tolSimilarity = 0.05;
	gaParam.similarityDebugLevel = 2;
	gaParam.energyReturnBfgs = -1.0e99;
	gaParam.removeSimilarStructures = 1;
	gaParam.bestIndividualSize = 3;

	// SYSTEM PARAMETERS
	int nAtoms = gaParam.numberOfParameters / 15;
	gaParam.nAtomTypes1 = nAtoms;
	gaParam.nAtomTypes2 = 0;
	gaParam.nAtomTypes3 = 0;
	gaParam.interactionPotentialType = 0;
	gaParam.generationToChangeInteraction = -1;
	for (int i = 0; i < nAtoms; i++)
		gaParam.atomLabels.push_back("H");

	// EXPERIMENT
	gaParam.isToMakeExperiment = false;
	gaParam.experimentFinalEnergy = 0;
	gaParam.experimentMaxSteps = 0;

	// INITIAL OPERATORS
	int nOperators = 5;
	double nOperatorsRate = 1.0e0 / (double)nOperators;
	gaParam.initialCreationRate.resize(nOperators);
	for(int i = 0; i < nOperators; i++)
	{
		gaParam.initialCreationRate[i] = nOperatorsRate;
	}


}




