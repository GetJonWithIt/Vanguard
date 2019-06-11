#ifndef MHDREDUCEDTESTS_H
#define MHDREDUCEDTESTS_H

#include "Solvers/mhdsecondordersolver.h"
#include <fstream>
using namespace std;

class MHDReducedTests
{
public:
    MHDReducedTests();

    static void solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solveDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void outputSolution(vector<MHDReducedStateVector> solution);
    static void outputSolution2D(vector<vector<MHDReducedStateVector> > solution);
};

#endif // MHDREDUCEDTESTS_H
