#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<iostream>

using namespace std;

int main()
{
	ifstream file_("creation-histogram.txt");
	int newIndividualsNumber = 10;
	int nOperators = 7;
	int nRuns = 50;
	int nMaxGenerations = 400;

	string line;
	vector< vector< vector<int> > > operatorsList;
	getline(file_,line);
	for(int i =0; i < nRuns; i++)
	{
		vector< vector<int> > runList;
		while(getline(file_,line))
		{
			stringstream convert;
			convert << line;
			string run;
			convert >> run;
			if(run == "run:")
				break;

			vector<int> operatorsLine(newIndividualsNumber);
			stringstream convert2;
			convert2 << line;
			for(int j = 0; j < operatorsLine.size(); j++)
				convert2 >> operatorsLine[j];
			
			
			vector<int> generationList(nOperators);
			for(int l = 0; l < generationList.size(); l++)
				generationList[l] = 0;


			for(int k = 0; k < operatorsLine.size(); k++)
				generationList[operatorsLine[k]]++;

			runList.push_back(generationList);

		}
		operatorsList.push_back(runList);
	}

	/* PRINT ONE
	for(int i = 0; i < operatorsList[0].size(); i++)
	{
		for(int j = 0; j < operatorsList[0][i].size();j++)
		{
			cout << operatorsList[0][i][j] << "  ";
		}
		cout << endl;
	}
	*/

	vector< vector<double> > meanRunList;
	int k = 0;
	while(k < nMaxGenerations)
	{
		double count = 0.0e0;
		vector<double> meanGenerationList(nOperators);
		for(int l = 0; l < nOperators; l++)
			meanGenerationList[l] = 0.0e0;

//		cout << "number fo runs: " << operatorsList.size() << endl;
		for(int i = 0; i < operatorsList.size(); i++)
		{
//			cout << "number of generation of run I   " << operatorsList[i].size() << endl;
//			cout << "k::  " << k << endl;
			if(operatorsList[i].size() > k)
			{
				for(int j = 0; j < nOperators; j++)
					meanGenerationList[j] += (double)operatorsList[i][k][j];
				count += 1.0e0;
			}
		}
		if(count == 0.0e0)
		{
			k++;
			continue;
		}

		for(int l = 0; l < meanGenerationList.size(); l++)
			meanGenerationList[l] /= count;

		meanRunList.push_back(meanGenerationList);
		k++;
	}

	ofstream printGraph_("meanRun.txt");
	for(int i = 0; i < meanRunList.size(); i++)
	{
		printGraph_ << i << "  ";
		for(int j = 0; j < meanRunList[i].size(); j++)
		{
			printGraph_ << meanRunList[i][j] << "  ";
		}
		printGraph_ << endl;
	}

	return 0;
}


