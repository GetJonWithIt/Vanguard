#ifndef MHDINTERMEDIATETESTS_H
#define MHDINTERMEDIATETESTS_H

#include "Solvers/mhdsecondordersolver.h"
#include <fstream>
using namespace std;

class MHDIntermediateTests
{
public:
    MHDIntermediateTests();

    static void solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solveDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void outputSolution(vector<MHDIntermediateStateVector> solution);
    static void outputSolution2D(vector<vector<MHDIntermediateStateVector> > solution);
};

#endif // MHDINTERMEDIATETESTS_H
