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

bool Similarity::checkSimilarity(std::vector<double> &x)
{
	if (method == 0)
	{
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
				if (printLevel > 1)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return true;
			}
		}
		return false;
	}
	else if (method == 1)
	{
		// Marques rmsd similarity
		vector<CoordXYZ> mol1(nAtoms);
		for (int i = 0; i < nAtoms; i++)
		{
			mol1[i].atomlabel = atomLabels[i];
			mol1[i].x = x[i];
			mol1[i].y = x[i + nAtoms];
			mol1[i].z = x[i + 2 * nAtoms];
		}
		MarquesEnantiomers mrq_;
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
			double distanceDiffererence = mrq_.marquesRmsd(mol1, mol2);
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 1)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return true;
			}
		}
		return false;
	}
	else if(method == -1)
	{
		return false;
	}
	else
	{
		cout << "SIMILARITY METHOD NOT FOUND" << endl;
		exit(1);
	}
}

//selected structures
bool Similarity::checkSimilarity(
	vector<double> &x,
	vector< vector<double> > &targedIndividuals
)
{
	if (method == 0)
	{
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
				return true;
			}
		}
		if (printLevel > 1)
		{
			if (targedIndividuals.size() != 0)
				pPrinting_->printSimilarityDistancesSelected(colectAllDifferences);
		}
		return false;
	}
	else if (method == 1)
	{
		vector<CoordXYZ> mol1(nAtoms);
		for (int i = 0; i < nAtoms; i++)
		{
			mol1[i].atomlabel = atomLabels[i];
			mol1[i].x = x[i];
			mol1[i].y = x[i + nAtoms];
			mol1[i].z = x[i + 2 * nAtoms];
		}

		vector<double> colectAllDifferences;
		MarquesEnantiomers mrq_;
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
			double distanceDiffererence = mrq_.marquesRmsd(mol1,mol2);
			colectAllDifferences.push_back(distanceDiffererence);
			if (distanceDiffererence < tolSimilarity)
			{
				if (printLevel > 1)
				{
					pPrinting_->printCreationIsEqual(
						distanceDiffererence,
						i);
				}
				return true;
			}
		}
		if (printLevel > 1)
		{
			if (targedIndividuals.size() != 0)
				pPrinting_->printSimilarityDistancesSelected(colectAllDifferences);
		}
		return false;
	}
	else if(method == -1)
        {
                return false;
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
		MarquesEnantiomers mrq_;
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
			double distanceDiffererence = mrq_.marquesRmsd(mol1, mol2);
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
		if (fitnessRank.size() != 0)
		{
			targetIndividuals.clear();
			for (int i = 0; i < bestIndividualsSize; i++)
				targetIndividuals.push_back(x_vec[fitnessRank[i]]);
		}
	}
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
	if ((printLevel > 0) &&
		(targetIndividuals.size() != 0) &&
		(activateIntoBfgsRmsd))
	{
		xToCheckBfgs.clear();
		pPrinting_->endlSimilarity();
	}
}

void Similarity::bestIndividualsCheck()
{
// checar se os melhores sao iguais
	if (targetIndividuals.size() != 0)
	{
		checkSimilarity(
			targetIndividuals[0],
			targetIndividuals);
	}

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

