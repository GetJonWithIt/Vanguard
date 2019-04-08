#ifndef MHDTESTS_H
#define MHDTESTS_H

#include "Solvers/mhdsecondordersolver.h"
#include <fstream>
using namespace std;

class MHDTests
{
public:
    MHDTests();

    static void solveDumbserTest1(int cellCount);
    static void solveDumbserTest2(int cellCount);
    static void solveDumbserTest3(int cellCount);

    static void outputSolution(vector<MHDStateVector> solution);
};

#endif // MHDTESTS_H
