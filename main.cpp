#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "BasicOperators.h"
#include "GeneticAlgorithm.h"
#include "Hooklaw.h"
#include "UserFitness.h"
#include "AuxMathGa.h"

#include "clusters/InitializeAtoms.h"
#include "clusters/ClustersFitness.h"
#include "clusters/ReadQuantumOutput.h"
#include "clusters/WriteQuantumInput.h"
#include "clusters/ReadGaInput.h"
#include "clusters/Fitness.h"
#include "clusters/ClustersOperators.h"
#include "clusters/Experiment.h"

using namespace std;
using namespace zerg;

/*
alteracoes do qga
- comando que pega um xyz, roda outra otimizacao e depois roda uma frequencia
!!!- se parar com o highlander, calcule a frequencia dele, se der bom, pare, se não SALVE esse cara e continue.
- portabilildade no fitness.cpp e no makefile (isnan)
- desativar leuituras desnecessarias - como frequencias

proximo paper
- algumas vezes eliminar não os piores, mas o highlander retirar os proximos a ele. o tempo e depois que o highlander comecar a se repetir muito.
- pensar sobre o highlander max iteration está 1.0e99.



estudo das formas de gerar clusters iniciais:
- Tem o de Sao Carlos, esfera, cubo e árvore.
- usar a algebra de polya para contar isomeros lennard jones dos bimetalicos,
  gerar todos e avaliar suas diferencas.
- colocar const no cluster operators por seguranca

- pop_size has to be multiple of four.

*/

void printAtomsVectorDouble(vector<double> & atoms, string testName = "teste.xyz");
void calculateMeanTestFormat(string name);
void generateExecutable(vector<string> argv);

int main(int argc, char *argv[])
{
	string experimentMethod;
	if (argc > 1)
	{
		stringstream convert0;
		convert0 << argv[1];
		convert0 >> experimentMethod;
	}
	vector<double> additionalParams;
	double aux1,aux2,aux3;
	stringstream convert1;
	int seed = 3;

	if(experimentMethod == "TwistOperator")
	{
		convert1 << argv[2] << "  " << argv[3] << "  " << argv[4];
		convert1 >> seed >> aux1 >> aux2;
		additionalParams.push_back(aux1);
		additionalParams.push_back(aux2);
		Experiment exp_;
		exp_.makeExperiment(seed, experimentMethod, additionalParams);
	}
	else if(experimentMethod == "GeometricCenterDisplacement")
	{
		convert1 << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5];
		convert1 >> seed >> aux1 >> aux2 >> aux3;
		additionalParams.push_back(aux1);
		additionalParams.push_back(aux2);
		additionalParams.push_back(aux3);
		Experiment exp_;
		exp_.makeExperiment(seed, experimentMethod, additionalParams);
	}
	else if(experimentMethod == "AutoAdjust")
	{
		convert1 << argv[2] << "  " << argv[3] << "  " << argv[4];
		convert1 >> seed >> aux1 >> aux2;
		additionalParams.push_back(aux1);
		additionalParams.push_back(aux2);
		Experiment exp_;
		exp_.makeExperiment(seed, experimentMethod, additionalParams);
	}
	else if (experimentMethod == "CalculateMean")
	{
		int size = 7;
		vector<string> name(size);
		for (size_t i = 0; i < size; i++)
		{
			stringstream convert;
			convert << (i + 1);

			name[i] = "cluster-result-"
				+ convert.str()
				+ ".csv";
		}
		//for (int i = 0; i < name.size(); i++)
		//calculateMeanTestFormat(name[i]);
		calculateMeanTestFormat("cluster-result-previousPaper.csv");
	}
	else if (experimentMethod == "experiment")
	{
		cout << "remove this histogram flag" << endl;
		exit(1);
		ofstream histogram_;
		histogram_.open("creation-histogram.txt", std::ofstream::out | std::ofstream::app);
		histogram_ << "run:  " << experimentMethod << "  seed:  " << seed << endl;
		histogram_.close();
		Experiment exp_;
		exp_.makeExperiment(seed, experimentMethod, additionalParams);
	}
	else
	{
		ReadGaInput readGa_;
		string gaInput;
		if(argc == 3)
			gaInput = argv[2];
		else
			gaInput = "GaInput.txt";

		readGa_.inputName = gaInput;

		readGa_.readGaInput();

		AuxMathGa::set_seed(readGa_.getSeed());

		zerg::GaParameters gaParam = readGa_.getGaParameters();

		vector<string> options = readGa_.getOptions();

		ClustersFitness clFit_(
			gaParam,
			options,
			readGa_.getGamessPath(),
			readGa_.getGamessScr(),
			readGa_.getGamessNprocess());

		GeneticAlgorithm ga1(clFit_, gaParam);

		ga1.ga_start();

		clFit_.printAllIndividuals("finalPopulation.xyz");
	}

	return 0;
}

void generateExecutable(vector<string> argv)
{
	stringstream conv;
	conv << argv[1] << "  " << argv[2] << "  " << argv[3];
	string auxName;
	double param1, param2;
	conv >> auxName >> param1 >> param2;
	string name = "./zerg.exe " + auxName + "  ";
	ofstream roda_("roda");
	roda_ << "#!/bin/bash" << endl;
	for (int i = 1; i <= 50; i++)
	{
		roda_ << name << i << "  " << param1
			<< "  " << param2 << endl;
	}
	roda_.close();
	system("chmod u+x roda");
}


void printAtomsVectorDouble(vector<double> & atoms, string testName)
{
	int natm = atoms.size() / 3;
	ofstream teste_(testName.c_str());
	teste_ << natm << endl << "t" << endl;
	for (int i = 0; i < natm; i++)
	{
		teste_ << "H "
			<< atoms[i] << "  "
			<< atoms[i + natm] << "  "
			<< atoms[i + 2 * natm] << endl;
	}
	teste_.close();
}

void calculateMeanTestFormat(string name)
{
	/*EXEMPLO
	int size = 3;
	vector<string> name(size);
	name[0] = "cluster-result-n-LargeAutoAdjust.csv";
	name[1] = "cluster-result-n-SmallAutoAdjust.csv";
	name[2] = "cluster-result-n-MediumAutoAdjust.csv";
	for (int i = 0; i < name.size(); i++)
		calculateMean(name[i]);
	*/

	ifstream read_(name.c_str());
	string auxline;
	vector<int> iterations;
	string a1, a2, a3, a4;
	int it;
	while (getline(read_, auxline))
	{
		if (auxline == "")
			break;
		stringstream line;
		line << auxline;
		line >> a1 >> a2 >> a3 >> a4 >> it;
		iterations.push_back(it);
	}
	read_.close();
	sort(iterations.begin(), iterations.end());
	string sName = name.erase(name.size() - 4, name.size());
	sName = sName.erase(0, 15);
	ofstream all_;
	all_.open("allMethods-Only.csv", std::ofstream::out | std::ofstream::app);
	if (iterations.size() == 0)
	{
		cout << "EMPTY FILE:  " << sName << endl;
		exit(1);
	}
	if (iterations[49] > 3000)
	{
		//all_ << sName << "  ;  " << "DISCARDED" << endl;
		//return;
	}
	int mean = 0;
	for (size_t i = 0; i < 50; i++)
	{
		mean += iterations[i];
	}
	all_ << sName << "  ;  " << (double)mean / 50.0e0 << endl;
}



