#ifndef HPRREDUCEDTESTS_H
#define HPRREDUCEDTESTS_H

#include "Solvers/hprsecondordersolver.h"
#include <fstream>
using namespace std;

class HPRReducedTests
{
public:
    HPRReducedTests();

    static void solveBartonTest1(int cellCount, int reinitialisationFrequency);

    static void solveZhangTest(int cellCount, int reinitialisationFrequency);

    static void outputSolution(vector<HPRReducedStateVector> solution);
};

#endif // HPRREDUCEDTESTS_H
