#ifndef HPRTESTS_H
#define HPRTESTS_H

#include "Solvers/hprsecondordersolver.h"
#include <fstream>
using namespace std;

class HPRTests
{
public:
    HPRTests();

    static void solveStokesFirstProblem(int cellCount, int subcyclingIterations);
    static void solveHeatConductionProblem(int cellCount, int subcyclingIterations);

    static void solveBartonTest1(int cellCount, int subcyclingIterations);

    static void solveZhangElasticPlasticTest(int cellCount, int subcyclingIterations);

    static void solve2DBartonTest1(int cellCount, int subcyclingIterations);

    static void solve2DZhangElasticPlasticTest(int cellCount, int subcyclingIterations);

    static void outputSolution(vector<HPRStateVector> solution, HPRMaterialParameters materialParameters);
    static void outputSolution2D(vector<vector<HPRStateVector> > solution, HPRMaterialParameters materialParameters);
};

#endif // HPRTESTS_H
