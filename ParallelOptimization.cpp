#include "ParallelOptimization.h"

#include <vector>

#include "Population.h"

using namespace std;
using namespace zerg;

namespace zerg{
void ParallelOptimization::do_optimization(Population &pop, std::vector<int> &dead_individuals)
{
	int size = dead_individuals.size();

	for(int i=0; i<size; i++)
	{
		pop.local_optimization(dead_individuals[i]);
	}
}

}







/*
#include <thread>
using std::thread;
#include <stdlib.h>
#include <Windows.h>
#include <vector>
using std::vector;


void fala_3seg()
{
	Sleep(3000);
	cout << "Fala de 3 segundos" << endl;
}

void fala_5seg()
{
	Sleep(5000);
	cout << "Fala de 5 segundos" << endl;
}

void fala_7seg()
{
	Sleep(7000);
	cout << "Fala de 7 segundos" << endl;
}


int main()
{
	thread vect[3];
	vect[0] = thread(fala_3seg);
	vect[1] = thread(fala_5seg);
	vect[2] = thread(fala_7seg);

	vect[0].join();
	vect[0] = thread(fala_3seg);


	vect[0].join();

	vect[1].join();
	vect[2].join();

	cin.get();
}
*/