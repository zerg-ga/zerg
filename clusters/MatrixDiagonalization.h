#ifndef MATRIXDIAGONALIZATION_H
#define MATRIXDIAGONALIZATION_H

#include <vector>

class MatrixDiagonalization
{
public:
	MatrixDiagonalization(){};
	~MatrixDiagonalization();
	void diagonalization(std::vector< std::vector<double> > &entering_matrix, int methodIn, int nElectrons);
	void symmetrizeEntryMatrix(std::vector< std::vector<double> > &matriz_entrada);
	std::vector< std::vector<double> > getEigenvectors();

	inline double getEigenvalueI(int i){ return eigenvalues[i]; }

private:
	double sign(double entryNumber);
	void choose_diagonalization_method(std::vector< std::vector<double> > &entering_matrix, int method, int nElectrons);
	std::vector<int> vector_ordering(std::vector<double> &vector_to_be_ordered);
	void vector_ordering_with_instructions(std::vector< std::vector<double> > &entering_vector,const std::vector<int> &vector_of_instructions);
	std::vector< std::vector<double> > eigenvectors;
	std::vector<double> eigenvalues;

	//pseudo
	int method;
	int nElectrons;
	int size;
	void modifiedFock(const std::vector< std::vector<double> > fockMatrix, int nElectrons_in);
	void rotations(std::vector< std::vector<double> > &modifiedMatrix );
	std::vector< std::vector<double> > pseudoEigenvectors;
	void twoVectorCombination(int i, int j, double alpha, double beta, std::vector< std::vector<double> > &modifiedMatrix);
	void testPseudo(const std::vector< std::vector<double> > fockMatrix, int nElectrons_in);

	void print_matrix(const std::vector< std::vector<double> > &entryMatrix);
};

#endif