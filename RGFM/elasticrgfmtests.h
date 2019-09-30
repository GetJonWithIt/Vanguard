#ifndef ELASTICRGFMTESTS_H
#define ELASTICRGFMTESTS_H

#include "elasticrgfmsolver.h"
#include <fstream>
using namespace std;

class ElasticRGFMTests
{
public:
    ElasticRGFMTests();

    static void solveBartonTest1Tilde(int cellCount, int subcyclingIterations, int reinitialisationFrequency);
    static void solveBartonTest1Star(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void solveZhangTestTilde(int cellCount, int subcyclingIterations, int reinitialisationFrequency);
    static void solveZhangTestStar(int cellCount, int subcyclingIterations, int reinitialisationFrequency);

    static void outputSolution(ElasticMultimaterialSystem multimaterialSystem);
};

#endif // ELASTICRGFMTESTS_H
