#ifndef MHDTESTS_H
#define MHDTESTS_H

#include "Solvers/mhdsecondordersolver.h"
#include <fstream>
using namespace std;

class MHDTests
{
public:
    MHDTests();

    static void solveDumbserTest1(int cellCount, int subcyclingIterations);
    static void solveDumbserTest2(int cellCount, int subcyclingIterations);
    static void solveDumbserTest3(int cellCount, int subcyclingIterations);

    static void solve2DDumbserTest1(int cellCount, int subcyclingIterations);
    static void solve2DDumbserTest2(int cellCount, int subcyclingIterations);
    static void solve2DDumbserTest3(int cellCount, int subcyclingIterations);

    static void outputSolution(vector<MHDStateVector> solution);
    static void outputSolution2D(vector<vector<MHDStateVector> > solution);
};

#endif // MHDTESTS_H
