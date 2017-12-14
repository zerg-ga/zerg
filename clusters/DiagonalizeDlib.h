// Fazer uma diagonaliza generica que herda dessa classe.
// eigen e o nome da biblioteca.
#ifndef DIAGONALIZEDLIB_H
#define DIAGONALIZEDLIB_H

#include <vector>

class DiagonalizeDlib
{
public:
	DiagonalizeDlib(std::vector< std::vector<double> > &entryMatrix);
	~DiagonalizeDlib();

	std::vector< std::vector<double> > eigenvectorsDlib;
	std::vector<double> eigenvaluesDlib;

private:
	//void diagonalizeWithDlib(std::vector< std::vector<double> > &entryMatrix);

};

#endif