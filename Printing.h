#ifndef PRINTING_H
#define PRINTING_H

#include <vector>
#include <fstream>
#include <string>

#include "StructOptions.h"

// WARNING
// DYNAMIC ALLOCATION, DON'T USE EXCEPTION HANDLING

namespace zerg {
	class Printing
	{
	public:
		Printing();

		~Printing();


		//main output - Read input debug
		void showExperimentMethod(std::string experimentMethod);
		void gaStartInputReading();
		void showInputLines(std::string auxline);
        	void showAllParameters(
                	zerg::GaParameters gaParam,
                	std::string gamessNproc,
               		std::string projectName,
                	std::string interactionPotential,
                	std::string gamessPath,
                	std::string gamessScr,
                	std::string gamessHeader,
                	std::vector<std::string> baseFiles);
		void endOfGamessOptions();

		//main output - GeneticAlgorithm
		void writeOpenMessage();
		void generationMessage(int gen);
		void energyMessage(
			std::vector<int> &fitness_rank, 
			std::vector<int> &dead_individuals,
			std::vector<double> & fitness_energies);
		void generationEndMessage();
		void highlanderMessage(int i, double frequency);
		void lastIndividualMessage(double frequency);
		void endMessage();

		// Creation output
		void setCreationDebug(int creationDebug);
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
