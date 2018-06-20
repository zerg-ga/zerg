#include "Similarity.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "Coordstructs.h"
#include "../Printing.h"
#include "MarquesEnantiomers.h"

using namespace std;
using namespace zerg;

Similarity::Similarity()
{
	bestIndividualsSize = 3;
	printLevel = 2;
	activateIntoBfgsRmsd = true;
	energyReturnBfgs = -1.0e99;
	bfgsSteps = 0;
}

Similarity::~Similarity(){}

void Similarity::startSimilarity(
	zerg::GaParameters & gaParam,
	Printing * pPrinting_in,
	Random * rand_in)
{
	activateIntoBfgsRmsd = gaParam.activateIntoBfgs;
	printLevel = gaParam.similarityDebugLevel;
	method = gaParam.similarityMethod;
	nAtoms = gaParam.numberOfParameters / 3;
	energyReturnBfgs = gaParam.energyReturnBfgs;
	tolSimilarity = gaParam.tolSimilarity;
	maxDistance = gaParam.maxDistance;
	minDistance = gaParam.minDistance;
	pPrinting_ = pPrinting_in;
	atomLabels = gaParam.atomLabels;
	bestIndividualsSize = gaParam.bestIndividualSize;

	if (printLevel != 0)
		pPrinting_->openSimilarityFile();

	marqRmsd_.setSeed(rand_in);
}


bool Similarity::checkLimitations(std::vector<double> &x)
{
	tempDistance.clear();
	tempDistance = calcAndSortAllDistances(x);
	if (tempDistance[tempDistance.size() - 1] > maxDistance)
		return true;
	if (tempDistance[0] < minDistance)
		return true;
	return false;

/* If I need to compare all to all - use this
	if (method == 1)
	{
		xTemp.clear();
		xTemp = x;
	}
*/
}

int Similarity::checkSimilarity(std::vector<double> &x)
{
	if (method == 0)
	{
		if (printLevel > 1)
			pPrinting_->printSimInitOutBfgs();

		// Distance similarity
		int size = tempDistance.size();
		for (size_t i = 0; i < allDistances.size(); i++)
		{
			double distanceDiffererence = 0.e0;
			for (int k = 0; k < size; k++)
				distanceDiffererence += abs(tempDistance[k] - allDistances[i][k]);

			distanceDiffererence /= (double)size;
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 0)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return i;
			}
		}
		if (printLevel > 1)
			pPrinting_->printSimUnique();
		return -1;
	}
	else if (method == 1)
	{
		if (printLevel > 1)
			pPrinting_->printSimInitOutBfgs();

		// Marques rmsd similarity
		vector<CoordXYZ> mol1(nAtoms);
		for (int i = 0; i < nAtoms; i++)
		{
			mol1[i].atomlabel = atomLabels[i];
			mol1[i].x = x[i];
			mol1[i].y = x[i + nAtoms];
			mol1[i].z = x[i + 2 * nAtoms];
		}
		for (size_t i = 0; i < targetIndividuals.size(); i++)
		{
			vector<CoordXYZ> mol2(nAtoms);
			for (int k = 0; k < nAtoms; k++)
			{
				mol2[k].atomlabel = atomLabels[k];
				mol2[k].x = targetIndividuals[i][k];
				mol2[k].y = targetIndividuals[i][k + nAtoms];
				mol2[k].z = targetIndividuals[i][k + 2 * nAtoms];
			}
			double distanceDiffererence = marqRmsd_.marquesRmsd(mol1, mol2);
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 0)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return i;
			}
		}
		if (printLevel > 1)
			pPrinting_->printSimUnique();
		return -1;
	}
	else if(method == -1)
	{
		return -1;
	}
	else
	{
		cout << "SIMILARITY METHOD NOT FOUND" << endl;
		exit(1);
	}
}

//selected structures
int Similarity::checkSimilarity(
	vector<double> &x,
	vector< vector<double> > &targedIndividuals
)
{
	if (method == 0)
	{
		if (printLevel > 1)
			pPrinting_->printSimInitOutBfgs();

		// Distance similarity
		int size = tempDistance.size();
		vector<double> colectAllDifferences;
		for (size_t i = 0; i < targedIndividuals.size(); i++)
		{
			vector<double> distTargetI = calcAndSortAllDistances(targedIndividuals[i]);
			double distanceDiffererence = 0.e0;
			for (int k = 0; k < size; k++)
				distanceDiffererence += abs(tempDistance[k] - distTargetI[k]);

			distanceDiffererence /= (double)size;
			colectAllDifferences.push_back(distanceDiffererence);
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 1)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return i;
			}
		}
		if (printLevel > 0)
		{
			if (targedIndividuals.size() != 0)
				pPrinting_->printSimilarityDistancesSelected(colectAllDifferences);
		}
		if (printLevel > 1)
			pPrinting_->printSimUnique();
		return -1;
	}
	else if (method == 1)
	{
		if (printLevel > 1)
			pPrinting_->printSimInitOutBfgs();

		vector<CoordXYZ> mol1(nAtoms);
		for (int i = 0; i < nAtoms; i++)
		{
			mol1[i].atomlabel = atomLabels[i];
			mol1[i].x = x[i];
			mol1[i].y = x[i + nAtoms];
			mol1[i].z = x[i + 2 * nAtoms];
		}

		vector<double> colectAllDifferences;
		for (size_t i = 0; i < targedIndividuals.size(); i++)
		{
			vector<CoordXYZ> mol2(nAtoms);
			for (int k = 0; k < nAtoms; k++)
			{
				mol2[k].atomlabel = atomLabels[k];
				mol2[k].x = targedIndividuals[i][k];
				mol2[k].y = targedIndividuals[i][k + nAtoms];
				mol2[k].z = targedIndividuals[i][k + 2 * nAtoms];
			}
			double distanceDiffererence = marqRmsd_.marquesRmsd(mol1,mol2);
			colectAllDifferences.push_back(distanceDiffererence);
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 1)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return i;
			}
		}
		if (printLevel > 0)
		{
			if (targedIndividuals.size() != 0)
				pPrinting_->printSimilarityDistancesSelected(colectAllDifferences);
		}
		if (printLevel > 1)
			pPrinting_->printSimUnique();
		return -1;
	}
	else if(method == -1)
        {
                return -1;
        }
	else
	{
		cout << "SIMILARITY METHOD NOT FOUND" << endl;
		exit(1);
	}
}

void Similarity::saveXToCheckBfgs(std::vector<double> & x)
{
	xToCheckBfgs = x;
}

double Similarity::checkSimilarityIntoBfgs()
{
	if (!activateIntoBfgsRmsd)
		return 0.0e0;
	bfgsSteps++;
	vector<double> colectAllDifferences;
	if (method == 0)
	{
		// CONDICAO DE SIMILARIDADE DE DISTANCIAS
		int size = tempDistance.size();
		for (size_t i = 0; i < targetIndividuals.size(); i++)
		{
			vector<double> distTargetI = calcAndSortAllDistances(targetIndividuals[i]);
			double distanceDiffererence = 0.e0;
			for (int k = 0; k < size; k++)
				distanceDiffererence += abs(tempDistance[k] - distTargetI[k]);

			distanceDiffererence /= (double)size;
			colectAllDifferences.push_back(distanceDiffererence);
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 1)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return energyReturnBfgs;
			}
		}
		if (printLevel > 0)
		{
			if (targetIndividuals.size() != 0)
				pPrinting_->printSimilarityDistances(colectAllDifferences);
		}
		return 0.0e0;
	}
	else if (method == 1)
	{
		vector<CoordXYZ> mol1(nAtoms);
		for (int i = 0; i < nAtoms; i++)
		{
			mol1[i].atomlabel = atomLabels[i];
			mol1[i].x = xToCheckBfgs[i];
			mol1[i].y = xToCheckBfgs[i + nAtoms];
			mol1[i].z = xToCheckBfgs[i + 2 * nAtoms];
		}

		vector<double> colectAllDifferences;
		for (size_t i = 0; i < targetIndividuals.size(); i++)
		{
			vector<CoordXYZ> mol2(nAtoms);
			for (int k = 0; k < nAtoms; k++)
			{
				mol2[k].atomlabel = atomLabels[k];
				mol2[k].x = targetIndividuals[i][k];
				mol2[k].y = targetIndividuals[i][k + nAtoms];
				mol2[k].z = targetIndividuals[i][k + 2 * nAtoms];
			}
			double distanceDiffererence = marqRmsd_.marquesRmsd(mol1, mol2);

			colectAllDifferences.push_back(distanceDiffererence);
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 1)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return energyReturnBfgs;
			}
		}
		if (printLevel > 0)
		{
			if (targetIndividuals.size() != 0)
				pPrinting_->printSimilarityDistances(colectAllDifferences);
		}
		return 0.0e0;
	}
	else if(method == -1)
		return false;
	else
	{
		cout << "SIMILARITY METHOD NOT FOUND" << endl;
		exit(1);
	}
}


void Similarity::addTargetIndividuals(
	std::vector< std::vector<double> > & x_vec,
	std::vector<int> & fitnessRank)
{
	if (method == 1)
	{
		if(distinctStruct.size() != 0)
		{
			int size;
			if(distinctStruct.size() < bestIndividualsSize)
				size = distinctStruct.size();
			else
				size = bestIndividualsSize;

			targetIndividuals.clear();
			for (int i = 0; i < size; i++)
				targetIndividuals.push_back(x_vec[distinctStruct[i]]);

		}
		else if (fitnessRank.size() != 0)
		{
			targetIndividuals.clear();
			for (int i = 0; i < bestIndividualsSize; i++)
				targetIndividuals.push_back(x_vec[fitnessRank[i]]);
		}
	}
}


void Similarity::bestIndividualsCorrections(
	std::vector<int> & corrections)
{
	if(method == 1)
	{
		for(int i = 0; i < bestIndividualsSize - 1; i++)
		{
			for(int j = i + 1; j < bestIndividualsSize; j++)
			{
				vector< vector<double> > indJ(1);
				indJ[0] = targetIndividuals[j];
				int equal = checkSimilarity(
						targetIndividuals[i],
						indJ);
				if(equal != -1)
					corrections.push_back(j);
			}
		}
	}
}


void Similarity::compareAllIndividuals(
	std::vector< std::vector<double> > & x_vec,
	std::vector<int> & fitnessRank)
{
	if(method == 1)
	{
		distinctStruct.clear();
		distinctStruct.push_back(fitnessRank[0]);
		if(distinctStruct.size() == (size_t)bestIndividualsSize)
		{
			pPrinting_->printDistinctStructures(distinctStruct);
			return;
		}

		for(int i = 1; i < fitnessRank.size() - 1; i++)
		{
			bool foundEqual = false;
			for(int j = 0; j < distinctStruct.size(); j++)//check if any distinct is equal to i
			{
				vector< vector<double> > indI(1);
				vector< vector<double> > indJ(1);
				indI[0] = x_vec[fitnessRank[i]];
				indJ[0] = x_vec[distinctStruct[j]];
				int equal = checkSimilarity(
						indI[0],
						indJ);
				if(equal != -1)//equal to some other structure
				{
					x_vec[fitnessRank[i]] = x_vec[distinctStruct[j]];
					foundEqual = true;
					break;
				}
			}
			if(!foundEqual)
			{
				distinctStruct.push_back(fitnessRank[i]);
				if(distinctStruct.size() == (size_t)bestIndividualsSize)
				{
					pPrinting_->printDistinctStructures(distinctStruct);
					return;
				}
			}
		}
	}

	pPrinting_->printDistinctStructures(distinctStruct);


}



void Similarity::appendTosimilarity()
{
	if (method == 0)
		allDistances.push_back(tempDistance);



/* If compare all to all is needed use this
else if (method == 1)
		allCoordinates.push_back(xTemp);
*/
}

void Similarity::printNewBfgsInd()
{
	if ((targetIndividuals.size() != 0) &&
		(activateIntoBfgsRmsd))
	{
		xToCheckBfgs.clear();
	}
	if ((activateIntoBfgsRmsd) && (printLevel > 1))
		pPrinting_->endlSimilarity();
}

void Similarity::printBfgsSteps()
{
	if (activateIntoBfgsRmsd)
		pPrinting_->printBfgsSteps(bfgsSteps);
}
















std::vector<double> Similarity::calcAndSortAllDistances(std::vector<double> &x)
{
	vector<double> auxDistances;
	for (int i = 0; i < (nAtoms - 1); i++)
	{
		for (int j = (i + 1); j < nAtoms; j++)
		{
			double dist = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + nAtoms] - x[j + nAtoms])*(x[i + nAtoms] - x[j + nAtoms]) +
				(x[i + 2 * nAtoms] - x[j + 2 * nAtoms])*(x[i + 2 * nAtoms] - x[j + 2 * nAtoms]));
			auxDistances.push_back(dist);
		}
	}
	sort(auxDistances.begin(), auxDistances.end());
	return auxDistances;
}

std::vector<double> Similarity::calcAndSortDistancesOverI(std::vector<double> &x, int i)
{
	vector<double> auxDistances;
	for (int j = 0; j < nAtoms; j++)
	{
		if (i == j)
			continue;
		double dist = sqrt(
			(x[i] - x[j])*(x[i] - x[j]) +
			(x[i + nAtoms] - x[j + nAtoms])*(x[i + nAtoms] - x[j + nAtoms]) +
			(x[i + 2 * nAtoms] - x[j + 2 * nAtoms])*(x[i + 2 * nAtoms] - x[j + 2 * nAtoms]));
		auxDistances.push_back(dist);
	}
	sort(auxDistances.begin(), auxDistances.end());
	return auxDistances;
}

std::vector<double> Similarity::calcDistanceToCenter(vector<double> &x)
{
	vector<double> auxDistances(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		auxDistances[i] = sqrt(
			x[i] * x[i] +
			x[i + nAtoms] * x[i + nAtoms] +
			x[i + 2 * nAtoms] * x[i + 2 * nAtoms]);
	}
	return auxDistances;
}

double Similarity::calcDistancesOverIAndGetMin(vector<double> &x, int i)
{
	vector<double> auxDistances;
	for (int j = 0; j < nAtoms; j++)
	{
		if (i == j)
			continue;
		double dist = sqrt(
			(x[i] - x[j])*(x[i] - x[j]) +
			(x[i + nAtoms] - x[j + nAtoms])*(x[i + nAtoms] - x[j + nAtoms]) +
			(x[i + 2 * nAtoms] - x[j + 2 * nAtoms])*(x[i + 2 * nAtoms] - x[j + 2 * nAtoms]));
		auxDistances.push_back(dist);
	}
	return *min_element(auxDistances.begin(), auxDistances.end());
}

void Similarity::printEndIntoBfgs()
{
	if((activateIntoBfgsRmsd) && (printLevel > 1))
		pPrinting_->printEndIntoBfgs();
}

void Similarity::printTitle()
{
	if(printLevel > 1)
		pPrinting_->printSimTitle();
}

void Similarity::printInitialTitle()
{
	if(printLevel > 0)
		pPrinting_->printInitialSimTitle();
}

