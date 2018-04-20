#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>

#include "BasicOperators.h"
#include "GeneticAlgorithm.h"
#include "Hooklaw.h"
#include "UserFitness.h"
#include "AuxMathGa.h"
#include "Printing.h"
#include "Random.h"

#include "clusters/Hungarian.h"
#include "clusters/MarquesEnantiomers.h"
#include "clusters/InitializeAtoms.h"
#include "clusters/ClustersFitness.h"
#include "clusters/ReadQuantumOutput.h"
#include "clusters/WriteQuantumInput.h"
#include "clusters/ReadGaInput.h"
#include "clusters/Fitness.h"
#include "clusters/ClustersOperators.h"
#include "clusters/Experiment.h"
#include "clusters/GamessIntoBfgs.h"

using namespace std;
using namespace zerg;

void printAtomsVectorDouble(vector<double> & atoms, string testName = "teste.xyz");
void calculateMeanTestFormat(string name);
void generateExecutable(vector<string> argv);
vector<double> readXyz(string xyzName);

std::vector< std::vector<double> > globalInd;

int main(int argc, char *argv[])
{

	/*
	MarquesEnantiomers marques_;
	marques_.calculateMarquesEnantiomers(
                "cluster-al9cu2-1",
                "cluster-al9cu2-2");
	*/

	/*
	vector<string> options(5);
	options[0] = "mopac";
	options[1] = "teste-mopac";
	options[2] = "PM7";
	options[3] = "";
	options[4] = "";
        WriteQuantumInput writeInp_(options);
        int nAtoms = 3;
        vector<CoordXYZ> mol(nAtoms);
	mol[0].atomlabel = "H";
	mol[1].atomlabel = "H";
	mol[2].atomlabel = "H";
	mol[0].x = 0.0e0;
	mol[1].x = 1.0e0;
	mol[2].x = 2.0e0;
	mol[0].y = 0.0e0;
	mol[1].y = 0.0e0;
	mol[2].y = 0.0e0;
	mol[0].z = 0.0e0;
	mol[1].z = 0.0e0;
	mol[1].z = 0.0e0;

        writeInp_.createInput(mol);

	return 0;

	ReadQuantumOutput readQ_("mopac");
	readQ_.readOutput("teste.out");
	vector<CoordXYZ> mol2 = readQ_.getCoordinates();
	for(size_t i = 0; i < mol.size(); i++)
	{
		cout << mol2[i].atomlabel << "  "
			<< mol2[i].x << endl;
	}
	cout << "energy:  " << readQ_.getEnergy();
	

	return 0;
	*/

	Printing printing_;
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
	if(experimentMethod == "")
		printing_.showExperimentMethod("GA run with input");
	else
		printing_.showExperimentMethod(experimentMethod);

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
		size_t size = 7;
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
	else if (experimentMethod == "test")
	{
		ReadQuantumOutput readQ_("gamess");
		readQ_.readOutput("na4.log");
		vector<CoordXYZ> mol = readQ_.getCoordinates();
		for(size_t i =0; i < mol.size(); i++)
		{
			cout << mol[i].atomlabel << "  " << mol[i].x << endl;
		}
		cout << "passou " << mol.size() << endl;
		return 0;	
	}
	else if (experimentMethod == "frequency")
	{
		Printing * pPrinting_ = &printing_;
		ReadGaInput readGa_(pPrinting_);
		string xyzName;
		string gaInput;
		gaInput = "GaInput.txt";
		convert1 << argv[2];
		convert1 >> xyzName;

		readGa_.inputName = gaInput;
		readGa_.readGaInput();
		zerg::GaParameters gaParam = readGa_.getGaParameters();
		vector<string> options = readGa_.getOptions();
		if (options.size() == 0)
		{
			cout << "Gamess options not found - check GaInput" << endl;
			exit(1);
		}
		Fitness fit_;
		vector<double> xCoordinates = readXyz(xyzName);
		options[1] = xyzName + "-";
		double energy = fit_.runGamess(
			xCoordinates,
			options,
			readGa_.getGamessPath(),
			readGa_.getGamessScr(),
			readGa_.getGamessNprocess());

		double frequency = fit_.runGamessFrequency(
			0,
			xCoordinates,
			options,
			readGa_.getGamessPath(),
			readGa_.getGamessScr(),
			readGa_.getGamessNprocess());

		cout << "energy: " << energy << "  frequency: " << frequency << endl;
	}
	else
	{
		ReadGaInput readGa_(&printing_);
		string gaInput;
		if(argc == 3)
			gaInput = argv[2];
		else
			gaInput = "GaInput.txt";

		readGa_.inputName = gaInput;

		readGa_.readGaInput();

		srand(readGa_.getSeed());

		Random rand_;

		zerg::GaParameters gaParam = readGa_.getGaParameters();

		vector<string> options = readGa_.getOptions();

		ClustersFitness clFit_(
			&rand_,
			gaParam,
			options,
			readGa_.getGamessPath(),
			readGa_.getGamessScr(),
			readGa_.getGamessNprocess(),
			&printing_);

		GeneticAlgorithm ga1(
			clFit_, 
			gaParam,
			&rand_,
			&printing_);

		ga1.ga_start();

		clFit_.printAllIndividuals("finalPopulation.xyz");

		clFit_.printBfgsSteps();


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

vector<double> readXyz(string xyzName)
{
	ifstream xyzFile_(xyzName.c_str());
	string line;
	stringstream convert;
	int natoms;
	getline(xyzFile_, line);
	convert << line;
	convert >> natoms;
	getline(xyzFile_, line);
	
	vector<double> atomX;
	vector<double> atomY;
	vector<double> atomZ;
	for (int i = 0; i < natoms; i++)
	{
		getline(xyzFile_, line);
		stringstream convert1;
		convert1 << line;
		string flag;
		double ax, ay, az;
		convert1 >> flag >> ax >> ay >> az;
		atomX.push_back(ax);
		atomY.push_back(ay);
		atomZ.push_back(az);
	}
	vector<double> coordinates(3 * natoms);
	for (int i = 0; i < natoms; i++)
	{
		coordinates[i] = atomX[i];
		coordinates[i + natoms] = atomY[i];
		coordinates[i + 2 * natoms] = atomZ[i];
	}
	return coordinates;
}



