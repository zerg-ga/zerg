#ifndef PRINTING_H
#define PRINTING_H

#include <vector>
#include <fstream>

// WARNING
// DYNAMIC ALLOCATION, DON'T USE EXCEPTION HANDLING

namespace zerg {
	class Printing
	{
	public:
		Printing(int creationDebug);

		~Printing();

		//main output - GeneticAlgorithm
		void writeOpenMessage();
		void generationMessage(int gen);
		void energyMessage(
			std::vector<int> &fitness_rank, 
			std::vector<int> &dead_individuals,
			std::vector<double> & fitness_energies);
		void generationEndMessage();
		void highlanderEndMessage();
		void endMessage();

		// Creation output
		void similarityProblem(int method);
		void setNewIndividualsError();
		void variationOfEachMethod(int method, double variation);
		void allVariations(std::vector<double> & methodMean);
		void factorToIncreaseDecrease();
		void printFactor(int method, double factor);
		void normalizedCreationRate(std::vector<double> & creationRate);
		void histogramTitle(int seed);
		void histogramPrint(int method);
		void histogramEndl();



	private:
		std::ofstream mainOutput_;

		int creationOutputDebugLevel; // 0-none ; 1-graph ; 2-all
		std::ofstream creationOutput_;

	};

}

#endif
