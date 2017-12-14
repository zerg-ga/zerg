#include "MatrixDiagonalization.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <stdlib.h>

#include "DiagonalizeDlib.h"

using namespace std;

MatrixDiagonalization::~MatrixDiagonalization(){}

void MatrixDiagonalization::diagonalization(std::vector< vector<double> > &enteringMatrix, int method, int nElectrons)
{
	symmetrizeEntryMatrix(enteringMatrix);
	choose_diagonalization_method(enteringMatrix, method, nElectrons);
}

void MatrixDiagonalization::choose_diagonalization_method(vector< vector<double> > &matriz_entrada, int methodIn, int nElectrons)
{
	method = methodIn;
	if (method == 0)
	{
		DiagonalizeDlib diagDlib_(matriz_entrada);
		eigenvalues = diagDlib_.eigenvaluesDlib;
		eigenvectors = diagDlib_.eigenvectorsDlib;
	}
	else if (method == 2)
	{
//		DiagonalizaEigen diag_eigen(matriz_entrada);
//		eigenvalues = diag_eigen.autovalores_eigen;
//		eigenvectors = diag_eigen.autovetores_eigen;
	}
	else if (method == 1)
	{
		modifiedFock(matriz_entrada, nElectrons);
		//pseudo diagonalization
	}
	else
	{
		cout << "ERRO NA ESCOLHA DO METODO DE DIAGONALIZAO\n"
			<< "OPCAO INDISPONIVEL\n";
		exit(2);
	}

	vector<int> instrucao_organiza = vector_ordering(eigenvalues);
	vector_ordering_with_instructions(eigenvectors, instrucao_organiza);

}

vector<int> MatrixDiagonalization::vector_ordering(vector<double> &vetor_entrada)
{
	int tamanho_vetor = vetor_entrada.size();
	vector<int> operacoes_organizam;

	double aux_menor;
	int aux_pos_menor;
	for (int i = 0; i < (tamanho_vetor - 1); i++)
	{

		aux_menor = vetor_entrada[i];
		aux_pos_menor = i;
		for (int j = i + 1; j < tamanho_vetor; j++)
		{
			if (aux_menor > vetor_entrada[j])
			{
				aux_menor = vetor_entrada[j];
				aux_pos_menor = j;
			}
		}

		// Se nao, nem precisa anotar a troca.
		if (aux_pos_menor != i)
		{
			vetor_entrada[aux_pos_menor] = vetor_entrada[i];
			vetor_entrada[i] = aux_menor;
			operacoes_organizam.push_back(i);
			operacoes_organizam.push_back(aux_pos_menor);
		}

	}

	// Nessa saida esta gravado todas as informacoes que foram feitas para organizar
	// o vetor;
	return operacoes_organizam;
}

// Rotina que organiza os autovetores em ordem crescente de autovalores.
void MatrixDiagonalization::vector_ordering_with_instructions(vector< vector<double> > &vetor_entrada, const vector<int> &vetor_organiza)
{
	int tamanho_vetor_organiza = vetor_organiza.size();
	int tamanho_vetor_entrada = vetor_entrada.size();
	int aux_pos1, aux_pos2;
	double aux_primeiro;

	for (int i = 0; i<tamanho_vetor_organiza; i += 2)
	{
		// o primeiro e igual ao segundo e o segundo e igual ao primeiro.
		aux_pos1 = vetor_organiza[i];
		aux_pos2 = vetor_organiza[i + 1];
		for (int j = 0; j < tamanho_vetor_entrada; j++)
		{
			aux_primeiro = vetor_entrada[j][aux_pos1];
			vetor_entrada[j][aux_pos1] = vetor_entrada[j][aux_pos2];
			vetor_entrada[j][aux_pos2] = aux_primeiro;
		}
	}

}


void MatrixDiagonalization::symmetrizeEntryMatrix(vector< vector<double> > &matriz_entrada)
{
	int n_bases = matriz_entrada.size();
	for (int i = 0; i < n_bases; i++)
	{
		for (int j = 0; j < n_bases; j++)
		{
			if (j < i)
			{
				matriz_entrada[i][j] = matriz_entrada[j][i];
			}
		}
	}
}


void MatrixDiagonalization::modifiedFock(const vector< vector<double> > fockMatrix, int nElectrons_in)
{
	nElectrons = nElectrons_in;
	double auxsoma;
	size = fockMatrix.size();
	int lumo = (size-nElectrons);// vector position

	vector< vector<double> > modifiedMatrixTemp;// nelec linhas e size colunas
	vector< vector<double> > modifiedMatrix; //nelec linhas e lumo colunas
	modifiedMatrixTemp.resize(nElectrons);
	modifiedMatrix.resize(nElectrons);
	for (int ii = 0; ii < (nElectrons); ii++)
	{
		modifiedMatrixTemp[ii].resize(size);
		modifiedMatrix[ii].resize(size-nElectrons);
	}
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < lumo; j++)
		{
			auxsoma = 0.0e0;
			for (int k = 0; k < size; k++)
			{
				auxsoma += eigenvectors[k][j] * fockMatrix[k][i];
			}
			modifiedMatrixTemp[j][i] = auxsoma;
		}
	}

	for (int i = 0; i < nElectrons; i++)
	{
		for (int j = lumo; j < size; j++)
		{
			auxsoma = 0.0e0;
			for (int k = 0; k < size; k++)
			{
				auxsoma += modifiedMatrixTemp[i][k] * eigenvectors[k][j];
			}
			modifiedMatrix[i][j-lumo] = auxsoma;

		}
	}
//	cout << "FIRST MODIFIED:  " << endl;
//	print_matrix(modifiedMatrix);

	rotations(modifiedMatrix);

//	testPseudo(fockMatrix, nElectrons);

}


// Routine just to test if pseudodiagonalization is correct
void MatrixDiagonalization::testPseudo(const vector< vector<double> > fockMatrix, int nElectrons_in)
{
	nElectrons = nElectrons_in;
	double auxsoma;
	size = fockMatrix.size();
	int lumo = (size - nElectrons);// vector position

	vector< vector<double> > modifiedMatrixTemp;// nelec linhas e size colunas
	vector< vector<double> > modifiedMatrix; //nelec linhas e lumo colunas
	modifiedMatrixTemp.resize(nElectrons);
	modifiedMatrix.resize(nElectrons);
	for (int ii = 0; ii < (nElectrons); ii++)
	{
		modifiedMatrixTemp[ii].resize(size);
		modifiedMatrix[ii].resize(size - nElectrons);
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < lumo; j++)
		{
			auxsoma = 0.0e0;
			for (int k = 0; k < size; k++)
			{
				auxsoma += pseudoEigenvectors[j][k] * fockMatrix[k][i];
			}
			modifiedMatrixTemp[j][i] = auxsoma;
		}
	}

	for (int i = 0; i < nElectrons; i++)
	{
		for (int j = lumo; j < size; j++)
		{
			auxsoma = 0.0e0;
			for (int k = 0; k < size; k++)
			{
				auxsoma += modifiedMatrixTemp[i][k] * pseudoEigenvectors[k][j];
			}
			modifiedMatrix[i][j - lumo] = auxsoma;

		}
	}

	vector< vector<double> > fockTemp1;// nelec linhas e size colunas
	vector< vector<double> > fockTemp2;// nelec linhas e size colunas
	fockTemp1 = fockMatrix;
	fockTemp2 = fockMatrix;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			fockTemp1[i][j] = 0.0e0;
			for (int k = 0; k < size; k++)
			{
				fockTemp1[i][j] += pseudoEigenvectors[k][i] * fockMatrix[k][j];
			}
		}
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			fockTemp2[i][j] = 0.0e0;
			for (int k = 0; k < size; k++)
			{
				fockTemp2[i][j] += fockTemp1[i][k] * pseudoEigenvectors[k][j];
			}
		}
	}


//		cout << "TEST MODIFIED:  " << endl;
//		print_matrix(modifiedMatrix);
//		cout << "DIAGONAL:  " << endl;
//		print_matrix(fockTemp2);
//		cout << "eigenvectors:  " << endl;
//		print_matrix(eigenvectors);
//		cout << "pseudo eigenvectors:  " << endl;
//		print_matrix(pseudoEigenvectors);

}



void MatrixDiagonalization::rotations(vector< vector<double> > &modifiedMatrix)
{
	int lumo = size - nElectrons;
	double tiny = 1.0e-7;
	pseudoEigenvectors = eigenvectors;

	double a,b,c,d,e,alpha,beta;
	for (int i = 0; i < nElectrons; i++)
	{
		for (int j = lumo; j < size; j++)
		{
			c = modifiedMatrix[i][j-lumo];
			if (abs(c) > tiny)
			{
				a = eigenvalues[i];
				b = eigenvalues[j];
				d = a - b;
				e = sign(d)*sqrt(4.0e0*c*c + d*d);
				alpha = sqrt(0.5e0*(1.0e0+d/e));
				beta = -sign(c)*sqrt(1.0e0 - alpha*alpha);
				twoVectorCombination(i, j, alpha, beta, modifiedMatrix);
			}
		}
	}

}

void MatrixDiagonalization::twoVectorCombination(int i, int j, double alpha, double beta, vector< vector<double> > &modifiedMatrix)
{
	vector<double> tempCi;
	vector<double> tempCj;
	tempCi.resize(size);
	tempCj.resize(size);

	for (int k = 0; k < size; k++)
	{
		tempCi[k] = alpha*pseudoEigenvectors[k][i] + beta*pseudoEigenvectors[k][j];
		tempCj[k] = alpha*pseudoEigenvectors[k][j] - beta*pseudoEigenvectors[k][i];
	}

	for (int ii = 0; ii < size; ii++)
	{
		pseudoEigenvectors[ii][i] = tempCi[ii];
		pseudoEigenvectors[ii][j] = tempCj[ii];
	}
}

vector< vector<double> > MatrixDiagonalization::getEigenvectors()
{
	if (method == 0)
	{
		return eigenvectors;
	}
	else
	{
		return pseudoEigenvectors;
	}

}


void MatrixDiagonalization::print_matrix(const vector< vector<double> > &entryMatrix)
{
	//cout.precision(6);
	int lines = entryMatrix.size();
	int cols = entryMatrix[0].size();
	for (int i = 0; i < lines; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << setiosflags(ios::fixed) << setprecision(6) << setw(12) <<
				entryMatrix[i][j] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}

double MatrixDiagonalization::sign(double x)
{
	if (x > 0) return 1.0e0;
	if (x < 0) return -1.0e0;
	return 0.0e0;
}
