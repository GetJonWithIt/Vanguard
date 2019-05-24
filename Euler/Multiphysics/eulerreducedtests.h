#ifndef EULERREDUCEDTESTS_H
#define EULERREDUCEDTESTS_H

#include "Solvers/secondordersolver.h"
#include <fstream>
using namespace std;

class EulerReducedTests
{
public:
    EulerReducedTests();

    static void solveToroTest1(int cellCount, int reinitialisationFrequency);

    static void solveFedkiwTest(int cellCount, int reinitialisationFrequency);

    static void solve2DToroTest1(int cellCount, int reinitialisationFrequency);

    static void solve2DFedkiwTest(int cellCount, int reinitialisationFrequency);

    static void outputSolution(vector<EulerReducedStateVector> solution);
    static void outputSolution2D(vector<vector<EulerReducedStateVector> > solution);
};

#endif // EULERREDUCEDTESTS_H
