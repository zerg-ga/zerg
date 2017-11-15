#ifndef GAPRINTING_H
#define GAPRINTING_H

#include <vector>
#include <fstream>

namespace zerg{
class GAPrinting
{
public:
	static void print_vector(const std::vector<double> &vec_in);
	static void print_vector(const std::vector<int> &vec_in);
	static void print_vector(const std::vector<double> &vec_in, std::ofstream &currentOut);
	static void print_vector(const std::vector<int> &vec_in, std::ofstream &currentOut);

};

}

#endif
