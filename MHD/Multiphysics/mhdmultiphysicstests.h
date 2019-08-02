#ifndef MHDMULTIPHYSICSTESTS_H
#define MHDMULTIPHYSICSTESTS_H

#include "Solvers/mhdsecondordersolver.h"
#include <fstream>
using namespace std;


class MHDMultiphysicsTests
{
public:
    MHDMultiphysicsTests();

    static void solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void outputSolution(vector<MHDMultiphysicsStateVector> solution);
};

#endif // MHDMULTIPHYSICSTESTS_H
