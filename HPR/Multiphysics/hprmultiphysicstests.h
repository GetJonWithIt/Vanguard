#ifndef HPRMULTIPHYSICSTESTS_H
#define HPRMULTIPHYSICSTESTS_H

#include "Solvers/hprsecondordersolver.h"
#include <fstream>
using namespace std;

class HPRMultiphysicsTests
{
public:
    HPRMultiphysicsTests();

    static void solveBartonTest1(int cellCount, int reinitialisationFrequency);

    static void solveZhangTest(int cellCount, int reinitialisationFrequency);

    static void outputSolution(vector<HPRMultiphysicsStateVector> solution);
};

#endif // HPRMULTIPHYSICSTESTS_H
