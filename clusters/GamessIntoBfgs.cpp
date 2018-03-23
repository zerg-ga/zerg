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
#include "Coordstructs.h"

using namespace std;

GamessIntoBfgs::GamessIntoBfgs(){}

GamessIntoBfgs::~GamessIntoBfgs(){}

void GamessIntoBfgs::runGamess(std::string inputName, string gamessExec) 
{
	pid_t child_pid;
	child_pid = fork();
	if (child_pid >= 0)
    	{
        	if (child_pid == 0) 
		{
			vector<string> commands(1);
			commands[0] = inputName;
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
			vector<CoordXYZ> mol;
			while(true)
			{

				// conferir com o numero de atomos

				sleep(1);
			       	readQ_.readOutput(inputName + "-.log");
				cout << "NSERCH:  " << readQ_.getNserch() << endl;
			        mol = readQ_.getCoordinates();
 			        for(size_t i = 0; i < mol.size(); i++)
       				{
          				     cout << mol[i].atomlabel << "  "
                   				  << mol[i].x << "  "
                     				  << mol[i].y << "  "
                       				  << mol[i].z << endl;
		        	} 		
				if(!checkGamess())
					break;

				if(readQ_.getNserch() > 10)
					killGamess();

			}
			cout << "done" << endl;
			exit(1);

			sleep(5);
			cout << "reading again" << endl;
			cout << "NSERCH:  " << readQ_.getNserch() << endl;
		        readQ_.readOutput(inputName + "-.log");
		        mol = readQ_.getCoordinates();
 		        for(size_t i = 0; i < mol.size(); i++)
       			{
          		     cout << mol[i].atomlabel << "  "
              			  << mol[i].x << "  "
                      		<< mol[i].y << "  "
                        	  << mol[i].z << endl;
			}
			sleep(10); 				
			cout << "reading again" << endl;
			cout << "NSERCH:  " << readQ_.getNserch() << endl;
		        readQ_.readOutput(inputName + "-.log");
		        mol = readQ_.getCoordinates();
 		        for(size_t i = 0; i < mol.size(); i++)
        		{
          		     cout << mol[i].atomlabel << "  "
              			  << mol[i].x << "  "
              			  << mol[i].y << "  "
                     		  << mol[i].z << endl;
		        } 				
			killGamess();
			// continuacao do codigo depois do fork
		}
	}
	else 
	{
	        cout << "Error on creating a parallel process" << endl
			<< "This shouldn't have happen, please, contact developers." << endl;
	        exit(0);
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

void GamessIntoBfgs::killGamess()
{
	system("ps ax | grep gamess.00.x > pidGamess.ga");
	ifstream readPid_("pidGamess.ga");
	string auxline;
	getline(readPid_, auxline);
	string pidGamess;
	readPid_ >> pidGamess;
	system(("kill " + pidGamess).c_str());
	readPid_.close();
	remove("pidGamess.ga");
	sleep(0.5);
}


