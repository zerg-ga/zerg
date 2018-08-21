#ifndef GAMESSINTOBFGS_H
#define GAMESSINTOBFGS_H

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <signal.h>
#include <string>
#include <vector>
#include <fstream>

#include "ReadQuantumOutput.h"
#include "Coordstructs.h"
#include "Similarity.h"

class GamessIntoBfgs
{
public:
	GamessIntoBfgs();
	~GamessIntoBfgs();

	double runGamess(
		int nAtoms,
		std::vector<CoordXYZ> & mol,
		std::string inputName, 
		std::string gamessExec,
		std::string nProc,
		Similarity * pSim_); 


private:
	bool checkGamess();
	
	void killGamess(pid_t child_pid);
	void killGamess();

};

#endif

