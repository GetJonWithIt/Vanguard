#ifndef EULERMULTIPHYSICSTESTS_H
#define EULERMULTIPHYSICSTESTS_H

#include "Solvers/secondordersolver.h"
#include <fstream>
using namespace std;

class EulerMultiphysicsTests
{
public:
    EulerMultiphysicsTests();

    static void solveToroTest1(int cellCount, int reinitialisationFrequency);

    static void solveFedkiwTest(int cellCount, int reinitialisationFrequency);

    static void solve2DToroTest1(int cellCount, int reinitialisationFrequency);

    static void solve2DFedkiwTest(int cellCount, int reinitialisationFrequency);

    static void outputSolution(vector<EulerMultiphysicsStateVector> solution);
    static void outputSolution2D(vector<vector<EulerMultiphysicsStateVector> > solution);
};

#endif // EULERMULTIPHYSICSTESTS_H
