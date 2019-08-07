#include "RunAtomistica.h"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdlib.h>

using namespace std;

vector<double> RunAtomistica::run(vector<double> & mol, double & energy)
{
        writeXyz(mol, "atomistica_run.xyz");
        system("python run_atomistica.py atomistica_run");
        vector<double> optimizedMol = readXyz("atomistica_run-done.xyz", energy);
        return optimizedMol;
}

void RunAtomistica::writeXyz(vector<double> & atoms, string testName)
{
        int natm = atoms.size() / 3;
        ofstream teste_;
        teste_.open(testName.c_str(), std::ofstream::out | std::ofstream::app);
        teste_ << natm << endl << "t" << endl;
        for (int i = 0; i < natm; i++)
        {
                teste_ << "C "
                        << atoms[i] << "  "
                        << atoms[i + natm] << "  "
                        << atoms[i + 2 * natm] << endl;
        }
        teste_.close();
}

vector<double> RunAtomistica::readXyz(string xyzName, double & energy)
{
        ifstream xyzFile_(xyzName.c_str());
        string line;
        stringstream convert;
        int natoms;
        getline(xyzFile_, line);
        convert << line;
        convert >> natoms;
        getline(xyzFile_, line);

        vector<double> mol(3 * natoms);
        for (int i = 0; i < natoms; i++)
        {
                getline(xyzFile_, line);
                stringstream convert1;
                convert1 << line;
                string label;
                double x,y,z;
                convert1 >> label
                        >> x
                        >> y
                        >> z;
                mol[i] = x;
                mol[i + natoms] = y;
                mol[i + 2*natoms] = z;
        }
        getline(xyzFile_, line);
        stringstream convert2;
        convert2 << line;
        convert2 >> energy;
        return mol;
}






