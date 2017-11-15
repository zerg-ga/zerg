#include <fstream>
#include <sstream>
#include <string>


using namespace std;

int main()
{
	string name = "./zerg.exe LargeAutoAdjust ";
	ofstream roda_("roda");
	roda_ << "#!/bin/bash" << endl;
	for(int i = 1; i <= 50; i++)
	{
		roda_ << name << i << endl;
	}
	return 0;


}


