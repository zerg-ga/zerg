#ifndef WRITEQUANTUMINPUT_H
#define WRITEQUANTUMINPUT_H

#include <fstream>
#include <string>
#include <vector>

#include "Coordstructs.h"

class WriteQuantumInput
{
public:
	WriteQuantumInput(std::vector<std::string> options);
	~WriteQuantumInput();

	std::string createInput(
		std::vector<CoordXYZ> &coordinates,
		int indexAddedAtFinalOfInputToSaveItsName = 0);

	// GLOBAL PUBLIC FUNCTIONS
	void changeProjectName(std::string newProjectName) { projectName = newProjectName; }

	//MOPAC PUBLIC FUNCTIONS
	void changeMopacHeader(std::string newHeader) { mopacHeader = newHeader; }

	//GAMESS PUBLIC FUNCIONS
	std::string getAtomName(int i);

private:
	//GLOBAL
	std::string type;
	void setInputProperties(std::vector<std::string> &options);
	std::string projectName;

	//MOPAC
	std::string mopacHeader;
	std::string externalParams;
	std::string centralMetal;
	void buildMopacInput(std::vector<CoordXYZ> &allAtoms, std::string inputName);

	//GAMESS
	std::vector<std::string> gamessHeader;
	std::vector<std::string> gamessAtomBasisFiles;
	std::vector<std::string> gamessEcpFiles;
	void readGamessAuxFiles();
	std::vector< std::string > atomName;

	std::vector< std::vector<std::string> > atomBasis;
	std::vector< std::vector<std::string> > atomEcp;
	void buildGamessInput(std::vector<CoordXYZ> &coordinates, std::string inputName);
	int getChargeFromAtom(std::string atomName);

};

#endif


/*
vector<string> options(16);
options[0] = "gamess";
options[1] = "teste";
options[2] = " $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE EXETYP=RUN MPLEVL=2 MAXIT=200 MULT=1";
options[3] = "   ISPHER=1 COORD=UNIQUE NOSYM=1 UNITS=ANGS $END";
options[4] = " $GUESS GUESS=HUCKEL $END";
options[5] = " $SYSTEM MWORDS=40 MEMDDI=20  $END";
options[6] = " $SCF DIRSCF=.FALSE. $END";
options[7] = " $DATA";
options[8] = "titulo";
options[9] = "C1";
options[10] = "EndOfHeader";
options[11] = "li-base.txt"; --- sempre dentro do auxFiles (end no final)
options[12] = "li-base.txt";
options[13] = "li-base.txt";
options[14] = "EndOfBasis";
options[15] = "NoEcp";
WriteQuantumInput writeInp_(options);
vector<CoordXYZ> mol(3);
mol[0].atomlabel = "Li";
mol[1].atomlabel = "Li";
mol[2].atomlabel = "Li";
mol[0].x = 0;
mol[0].y = 0;
mol[0].z = 0;
mol[1].x = 1;
mol[1].y = 0;
mol[1].z = 0;
mol[2].x = 2;
mol[2].y = 0;
mol[2].z = 0;
writeInp_.createInput(mol,5);


*/





