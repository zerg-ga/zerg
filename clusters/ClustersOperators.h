#ifndef CLUSTERSOPERATORS_H
#define CLUSTERSOPERATORS_H

#include <vector>
#include <string>

#include "InitializeAtoms.h"
#include "../BasicOperators.h"
#include "../StructOptions.h"

class ClustersOperators : public zerg::BasicOperators
{
public:
	ClustersOperators(int pop_size, int number_parameters);
	~ClustersOperators();

	//reimplement or use basic
	void startClustersOperators(zerg::GaParameters & gaParam);//final adjustments  (like crossover probability)

	bool create_individual(int creation_type, int target, int parent1, int parent2);
	bool operatorAdministration(int method, const std::vector<double> &operatorPerformance);// modify operators with it's performance
	bool check_similarity(int target);

	//keep
	virtual void local_optimization(int ind_i) = 0;
	virtual bool checkMinimum(int ind_i) { return true; }
	
protected:
	void appendTosimilarity(int ind_i);
	void translateToGeometricCenter(std::vector<double> & x);


	// OPERATORS PARAMETERS
	double scdo;
	double alfaMinGcdo;
	double alfaMaxGcdo;
	double wGcdo;
	double tetaMinTwisto;
	double tetaMaxTwisto;
	double contractionMinMtco;
	double contractionMaxMtco;


private:
	//set on startUserOperators()
	int nAtoms;
	double gamma;
	double rca;
	double adminLargeEnergyVariation;
	double maxDistance;

	bool sphereCutAndSplice(int target, int parent1, int parent2);

	std::vector<double> calculateRadius(std::vector<double> &x);

	std::vector<double> calcAndSortAllDistances(std::vector<double> &x);

	std::vector<double> calcAndSortDistancesOverI(std::vector<double> &x, int i);

	std::vector<double> calcDistanceToCenter(std::vector<double> &x);

	double calcDistancesOverIAndGetMin(std::vector<double> &x, int i);
	
	std::vector< std::vector<double> > allDistances;

	InitializeAtoms init_;

	void printAtomsVectorDouble(std::vector<double> & atoms, std::string testName);

	void translate(std::vector<double> & x, std::vector<double> & translateVector);

	//CLUSTERS OPERATORS
	std::vector<double> rondinaCartesianDisplacementOperator(const std::vector<double> & x);
	std::vector<double> rondinaGeometricCenterDisplacementOperator(const std::vector<double> & x);
	std::vector<double> rondinaTwistOperator(const std::vector<double> & x);
	std::vector<double> rondinaAngularOperator(const std::vector<double> & x);
	std::vector<double> rondinaAngularSurfaceOperator(const std::vector<double> & x);
	std::vector<double> rondinaMoveToCenterOperator(const std::vector<double> & x);
	std::vector<double> fredAngularSurfaceOperator(const std::vector<double> & x);
	std::vector<double> deavenHoCutSplice(const std::vector<double> & x1_parent, const std::vector<double> & x2_parent);


};

#endif

