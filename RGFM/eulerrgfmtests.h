#ifndef EULERRGFMTESTS_H
#define EULERRGFMTESTS_H

#include "rgfmsolver.h"
#include <fstream>
using namespace std;

class EulerRGFMTests
{
public:
    EulerRGFMTests();

    static void solveToroTest1Exact(int cellCount, int reinitialisationFrequency);
    static void solveToroTest1HLLC(int cellCount, int reinitialisationFrequency);

    static void solveFedkiwTestExact(int cellCount, int reinitialisationFrequency);
    static void solveFedkiwTestHLLC(int cellCount, int reinitialisationFrequency);

    static void solve2DToroTest1Exact(int cellCount, int reinitialisationFrequency);
    static void solve2DToroTest1HLLC(int cellCount, int reinitialisationFrequency);

    static void solve2DFedkiwTestExact(int cellCount, int reinitialistaionFrequency);
    static void solve2DFedkiwTestHLLC(int cellCount, int reinitialisationFrequency);

    static void outputSolution(MultimaterialSystem multimaterialSystem);
    static void outputSolution2D(MultimaterialSystem multimaterialSystem);
};

#endif // EULERRGFMTESTS_H
