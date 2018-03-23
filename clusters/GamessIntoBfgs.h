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

class GamessIntoBfgs
{
public:
	GamessIntoBfgs();
	~GamessIntoBfgs();

	void runGamess(std::string inputName, std::string gamessExec); 


private:
	bool checkGamess();
	
	void killGamess();

};

#endif

