#ifndef READGAINPUT_H
#define READGAINPUT_H

#include "../StructOptions.h"
#include "../Printing.h"

#include <string>
#include <vector>
#include <sstream>

class ReadGaInput
{
public:
	ReadGaInput(zerg::Printing * pPrinting_in);

	~ReadGaInput();

	void readGaInput();

	//void setExperimentDefaults(int seed);

	void setInputInformations(zerg::GaParameters gaParam_in);

	std::vector<std::string> getOptions();

	int getSeed() {return gaParam.seed;}

	zerg::GaParameters getGaParameters() { return gaParam; }

	std::string getGamessPath() { return gamessPath; }

	std::string getGamessNprocess() { return gamessNproc; }

	std::string getGamessScr() { return gamessScr; }

	std::string inputName;

private:
	int seed;
	std::vector<std::string> options;
	zerg::GaParameters gaParam;
	std::string gamessPath;
	std::string gamessScr;
	std::string gamessNproc;
	void setDefaults();
	void userDefinedSet(int userMethod);
	void baseFilesOrdering(
		std::vector<std::string> & baseFiles,
		int & n1,
		int & n2,
		int & n3);
	bool convertInputToBool(std::stringstream &convert);

	zerg::Printing * pPrinting_;

};


#endif

