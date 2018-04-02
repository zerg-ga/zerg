#define useDlib

#include "Fitness.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <stdlib.h>

#include "WriteQuantumInput.h"
#include "ReadQuantumOutput.h"
#include "Similarity.h"
#include "GamessIntoBfgs.h"

#ifdef useDlib
#include <dlib/optimization.h>
#include "FunctionDlib.h"
#include "FunctionDlibSim.h"
#include "DerivativeDlib.h"
#endif

using namespace std;

Fitness::Fitness()
{
	lowestPossibleEnergy = -1.0e90;
}

Fitness::~Fitness(){}

double Fitness::fit(vector<double> &point, int type)
{
	switch (type)
	{
	case 0:
		return lennardJones(point);
		break;

	case 1:
		cout << "ERROR ON: Fitness::fit - need to set parameters" << endl;
		exit(1);
		break;

	default:
		cout << "FITNESS FUNCTION NOT FOUND" << endl;
		exit(3);
	}
}

double Fitness::fit(
	vector<double> &point, 
	int type,
	const vector<double> params,
	const vector<int> atomTypes)
{
	switch (type)
	{
	case 0:
		return lennardJones(point);
		break;

	case 1:
		return gupta(point, params, atomTypes);
		break;

	default:
		cout << "FITNESS FUNCTION NOT FOUND" << endl;
		exit(3);
	}
}


double Fitness::gupta(
	vector<double> &x,
	const vector<double> atomsParameters,
	const vector<int> atomTypes)
{
	
	// x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3
	int natm = x.size() / 3;
	double r;
	double totalRep = 0.0e0;
	double totalAttr = 0.0e0;
	for (int i = 0; i < natm; i++)
	{
		double repSum = 0.0e0;
		double attrSum = 0.0e0;
		for (int j = 0; j < natm; j++)
		{
			if (j == i)
				continue;
			vector<double> parameters = getGuptaParameters(
				i,
				j,
				atomTypes,
				atomsParameters);

			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
				);

			// parameters: A, zeta, p, q e r0.

			r /= parameters[4];

			repSum += parameters[0]*exp(parameters[2]*(1.0e0 - r));

			attrSum += parameters[1]*parameters[1]*exp(2.0e0*parameters[3]*(1.0e0 - r));
		}
		totalRep += repSum;
		totalAttr += sqrt(attrSum);

	}

	return totalRep - totalAttr;

}

double Fitness::lennardJones(vector<double> &x)
{
	// x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3
	int natm = x.size() / 3;
	double r, r2, r4, r6, r12;
	double vlj = 0.0e0;
	for (int i = 0; i < (natm - 1); i++)
	{
		for (int j = (i + 1); j < natm; j++)
		{
			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
			);
			r2 = r * r;
			r4 = r2 * r2;
			r6 = r4 * r2;
			r12 = r6 * r6;
			vlj += 4.0e0 * (-1 / r6 + 1 / r12);
		}
	}
	//if (isnan(vlj))
	//	return 1.0e99;
	return vlj;
}



double Fitness::runGamess(
	vector<double> &x, 
	vector<string> &options, 
	string gamessPath,
	string gamessScr,
	string nProc)
{
	WriteQuantumInput writeInp_(options);
	int nAtoms = x.size() / 3;
	vector<CoordXYZ> mol(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		mol[i].atomlabel = writeInp_.getAtomName(i);
		mol[i].x = x[i];
		mol[i].y = x[i + nAtoms];
		mol[i].z = x[i + 2 * nAtoms];
	}

	writeInp_.createInput(mol);

	system(("rm " + gamessScr + "/" + options[1] + "*").c_str());

	system((gamessPath + "  " + options[1] + ".inp  00  " + nProc + " > " + options[1] + ".log").c_str());

	ReadQuantumOutput readQ_("gamess");

	readQ_.readOutput((options[1] + ".log").c_str());

	mol = readQ_.getCoordinates();
	if(mol.size() == 0)
	{
		return 0.0e0;
	}
	else
	{
		for (int i = 0; i < nAtoms; i++)
		{
			x[i] = mol[i].x;
			x[i + nAtoms] = mol[i].y;
			x[i + 2 * nAtoms] = mol[i].z;
		}
		return readQ_.getEnergy();
	}
}


double Fitness::runGamess(
	vector<double> &x, 
	vector<string> &options, 
	string gamessPath,
	string gamessScr,
	string nProc,
	Similarity * pSim_)
{
	WriteQuantumInput writeInp_(options);
	int nAtoms = x.size() / 3;
	vector<CoordXYZ> mol(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		mol[i].atomlabel = writeInp_.getAtomName(i);
		mol[i].x = x[i];
		mol[i].y = x[i + nAtoms];
		mol[i].z = x[i + 2 * nAtoms];
	}

	writeInp_.createInput(mol);

	system(("rm " + gamessScr + "/" + options[1] + "*").c_str());

	GamessIntoBfgs gameInto_;

	double energy = gameInto_.runGamess(
		mol.size(),
		mol,
		options[1] + ".inp",
		gamessPath,
		nProc,
		pSim_);

//	system((gamessPath + "  " + options[1] + ".inp  00  " + nProc + " > " + options[1] + ".log").c_str());
//	ReadQuantumOutput readQ_("gamess");
//	readQ_.readOutput((options[1] + ".inp-.log").c_str());
//	mol = readQ_.getCoordinates();
	if(mol.size() == 0)
	{
		return 0.0e0;
	}
	else
	{
		for (int i = 0; i < nAtoms; i++)
		{
			x[i] = mol[i].x;
			x[i + nAtoms] = mol[i].y;
			x[i + 2 * nAtoms] = mol[i].z;
		}
		return energy;
	}
}





double Fitness::runGamessFrequency(
	int numberOfOptimizations,
	vector<double> &x,
	vector<string> &options,
	string gamessPath,
	string gamessScr,
	string nProc)
{

	vector<string> frequencyOptions;

	//changing name
	stringstream convertInt;
	convertInt << numberOfOptimizations;
	string numberString;
	convertInt >> numberString;
	frequencyOptions.push_back(options[0]);
	string freqFileName = options[1] + "-highlander-" + numberString + "-frequency";
	frequencyOptions.push_back(freqFileName);

	string gamessHeader;
	for(size_t i = 0; i < options.size(); i++)
	{
		if(options[i] == "Header name")
		{
			gamessHeader = options[i+1];
			break;
		}
	}
	gamessHeader = "frequency-" + gamessHeader;

	ifstream gamHeader_(("auxFiles/" + gamessHeader).c_str());
	string auxline;
        while (getline(gamHeader_, auxline))
        {
		if (auxline.find("EndOfGamessHeader") != string::npos)
			break;
		frequencyOptions.push_back(auxline);
	}
	gamHeader_.close();
	frequencyOptions.push_back("EndOfHeader");

	int iEndHeader = -1;
	for(int i = 0; i < (int)options.size(); i++)
	{
		if(options[i] == "EndOfHeader")
		{
			iEndHeader = i;
		}
		if((i > iEndHeader) && (iEndHeader != -1))
		{
			frequencyOptions.push_back(options[i]);
		}
	}
		
	/*
	//replacing RUNTY=OPTIMIZE
	for (size_t i = 0; i < frequencyOptions.size(); i++)
	{
		if (frequencyOptions[i].find("RUNTYP=OPTIMIZE") != string::npos)
		{
			frequencyOptions[i].erase(frequencyOptions[i].find("RUNTYP=OPTIMIZE"), 16);
			frequencyOptions[i] += " RUNTYP=HESSIAN";
			break;
		}
	}
	*/

	WriteQuantumInput writeInp_(frequencyOptions);

	int nAtoms = x.size() / 3;
	vector<CoordXYZ> mol(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		mol[i].atomlabel = writeInp_.getAtomName(i);
		mol[i].x = x[i];
		mol[i].y = x[i + nAtoms];
		mol[i].z = x[i + 2 * nAtoms];
	}

	writeInp_.createInput(mol);
	
	system(("rm " + gamessScr + "/" + frequencyOptions[1] + "*").c_str());

	system((gamessPath + "  " + frequencyOptions[1] + ".inp  00  " + nProc + " > " + frequencyOptions[1] + ".log").c_str());

	ReadQuantumOutput readQ_("gamess");

	readQ_.activateDeactivateReadings("frequency", true);

	readQ_.activateDeactivateReadings("coordinates", false);

	readQ_.activateDeactivateReadings("energy", false);

	readQ_.readOutput((frequencyOptions[1] + ".log").c_str());

	return readQ_.getFirstFrequency();

}

double Fitness::optimizeEmpiricalPotential(
	std::vector<double> &x,
	int fitType,
	vector<double> &parameters,
	vector<int> &atomTypes,
	Similarity * pSim_)
{
#ifdef useDlib
	using namespace dlib;

	int size = x.size();
	column_vector starting_point(size);
	for (int i = 0; i < size; i++)
		starting_point(i) = x[i];

	double fMin = find_min(bfgs_search_strategy(),
		objective_delta_stop_strategy(1e-6),
		FunctionDlibSim(
			size, 
			fitType, 
			parameters,
			atomTypes,
			pSim_),
		DerivativeDlib(
			size, 
			fitType, 
			parameters,
			atomTypes),
		starting_point,
		lowestPossibleEnergy);

	if (fMin < lowestPossibleEnergy)
		fMin = 1.0e99;

	for (int i = 0; i < size; i++)
		x[i] = starting_point(i);

	return fMin;
#else
	return lennardJones(x);
#endif	
}



std::vector<double> Fitness::getGuptaParameters(
	int iPar,
	int jPar,
	const std::vector<int> atomTypes,
	const std::vector<double> atomsParameters)
{
	int maxElem = *max_element(atomTypes.begin(), atomTypes.end());

	vector<double> parameters(5);
	int k = -1;
	for (int i = 0; i <= maxElem; i++)
	{
		for (int j = i; j <= maxElem; j++)
		{
			k++;
			if (
				((atomTypes[iPar] == i) && (atomTypes[jPar] == j)) ||
				((atomTypes[iPar] == j) && (atomTypes[jPar] == i))
				)
			{
				for (int l = 0; l < 5; l++)
				{
					parameters[l] = atomsParameters[5 * k + l];
				}
				return parameters;
			}
		}
	}
	cout << "ERROR ON: Fitness::getGuptaParameters" << endl;
	exit(1);
}




double Fitness::runMopac(
	vector<double> &x, 
	vector<string> &options,
	vector<string> &atomLabels)
{
	WriteQuantumInput writeInp_(options);
	int nAtoms = x.size() / 3;
	vector<CoordXYZ> mol(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		mol[i].atomlabel = atomLabels[i];
		mol[i].x = x[i];
		mol[i].y = x[i + nAtoms];
		mol[i].z = x[i + 2 * nAtoms];
	}

	writeInp_.createInput(mol);

	system(("/opt/mopac/MOPAC2016.exe " + options[1] + ".mop ").c_str());

	ReadQuantumOutput readQ_("mopac");

	readQ_.readOutput((options[1] + ".out").c_str());

	mol = readQ_.getCoordinates();

	if(mol.size() == 0)
	{
		return 0.0e0;
	}
	else
	{
		for (int i = 0; i < nAtoms; i++)
		{
			x[i] = mol[i].x;
			x[i + nAtoms] = mol[i].y;
			x[i + 2 * nAtoms] = mol[i].z;
		}
		return readQ_.getEnergy();
	}
}

















/* EXMPLO DE OPTIMIZE
//main:
//InitializeAtoms init_;
//vector<double> x = init_.generateCluster(20, 0.2, 2.5);
//Fitness fit_;
//printAtomsVectorDouble(x, "teste1.xyz");
//fit_.optimizeEmpiricalPotential(x, 0);
//printAtomsVectorDouble(x, "teste2.xyz");
*/
