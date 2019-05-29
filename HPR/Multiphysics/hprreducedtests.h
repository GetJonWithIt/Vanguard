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

    static void solveZhangElasticPlasticTest(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DBartonTest1(int cellCount, int reinitialisationFrequency);

    static void solve2DZhangTest(int cellCount, int reinitialisationFrequency);

    static void solve2DZhangTest2(int cellCount, int reinitialisationFrequency);

    static void solve2DZhangElasticPlasticTest(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DZhangElasticPlasticTest2(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DHowellBallTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);
    static void solve2DHowellBallTest2(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void outputSolution(vector<HPRReducedStateVector> solution);
    static void outputSolution2D(vector<vector<HPRReducedStateVector> > solution);
};

#endif // HPRREDUCEDTESTS_H
