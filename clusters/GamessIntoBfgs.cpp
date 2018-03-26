#include "GamessIntoBfgs.h"

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
#include <sstream>

#include "ReadQuantumOutput.h"
#include "Similarity.h"
#include "Coordstructs.h"

using namespace std;

GamessIntoBfgs::GamessIntoBfgs(){}

GamessIntoBfgs::~GamessIntoBfgs(){}

double GamessIntoBfgs::runGamess(
	int nAtoms,
	std::vector<CoordXYZ> &mol,
	std::string inputName, 
	string gamessExec,
	string nProc,
	Similarity * pSim_) 
{
	pid_t child_pid;
	child_pid = fork();
	if (child_pid >= 0)
    	{
        	if (child_pid == 0) 
		{
			vector<string> commands(3);
			commands[0] = inputName;
			commands[1] = "00";
			commands[2] = nProc;
			remove((inputName + "-.log").c_str());
			int fd = open((inputName + "-.log").c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
			dup2(fd,1);
			dup2(fd,2);
			close(fd);
			const char *programname = gamessExec.c_str();
			const char **argv = new const char* [commands.size()+2];   
			argv[0] = programname;         
			for (size_t j = 0;  j < commands.size()+1;j++)
				argv[j+1] = commands[j].c_str();
			argv [commands.size()+1] = NULL;
				
			int status = execv(programname, (char **)argv);
			perror("execv");
			cout << "ERROR EXECUTING GAMESS:  " << status;
			exit(1);
		}
		else
		{
       			ReadQuantumOutput readQ_("gamess");
			int oldNserch = -1;
			mol.clear();
			while(true)
			{

				// conferir com o numero de atomos

				sleep(1);
			       	readQ_.readOutput(inputName + "-.log");
			        mol = readQ_.getCoordinates();
				if(((int)mol.size() == nAtoms) && (readQ_.getNserch() != oldNserch))
				{
					oldNserch = readQ_.getNserch();
					vector<double> x(3*mol.size());
					for(size_t i = 0; i < mol.size(); i++)
					{
						x[i] = mol[i].x;
						x[i + nAtoms] = mol[i].y;
						x[i + 2 * nAtoms] = mol[i].z;
					}
					pSim_->saveXToCheckBfgs(x);
					double simFlag = pSim_->checkSimilarityIntoBfgs();
					if(simFlag < -1.0e90)
						killGamess(child_pid);
		
					/*	
					cout << "NSERCH:  " << readQ_.getNserch() << endl;
			        	for(size_t i = 0; i < mol.size(); i++)
       					{
          					     cout << mol[i].atomlabel << "  "
                   					  << mol[i].x << "  "
                     					  << mol[i].y << "  "
                       					  << mol[i].z << endl;
		        		} 		
					*/

				}
				if(!checkGamess())
					break;
			}
			if(mol.size() == 0)
				return 0.0e0;
			else
				return readQ_.getEnergy();
		}
	}
	else 
	{
	        cout << "Error on creating parallel process" << endl
			<< "This shouldn't have happened, please, contact developers." << endl;
	        exit(1);
	}
}

bool GamessIntoBfgs::checkGamess()
{
	system("ps ax | grep gamess.00.x > pidGamess.ga");
	ifstream readPid_("pidGamess.ga");
	string auxline;
	getline(readPid_, auxline);
	readPid_.close();
	remove("pidGamess.ga");
	if(auxline.find("grep") != string::npos)
		return false;
	return true;
}

void GamessIntoBfgs::killGamess(pid_t child_pid)
{
	kill(child_pid, SIGTERM);
	system("ps ax | grep gamess.00.x > pidGamess.ga");
	ifstream readPid_("pidGamess.ga");
	string auxline;
	getline(readPid_, auxline);
	string pidGamess;
	readPid_ >> pidGamess;
	system(("kill " + pidGamess).c_str());
	readPid_.close();
	remove("pidGamess.ga");
	int status;
	wait(&status);
	sleep(0.5);
}


