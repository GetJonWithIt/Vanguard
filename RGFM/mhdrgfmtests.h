#ifndef MHDRGFMTESTS_H
#define MHDRGFMTESTS_H

#include "mhdrgfmsolver.h"
#include <fstream>
using namespace std;

class MHDRGFMTests
{
public:
    MHDRGFMTests();

    static void solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solveDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solve2DDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void outputSolution(MHDMultimaterialSystem multimaterialSystem);
    static void outputSolution2D(MHDMultimaterialSystem multimaterialSystem);
};

#endif // MHDRGFMTESTS_H
