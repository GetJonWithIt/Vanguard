#ifndef ELASTICREDUCEDTESTS_H
#define ELASTICREDUCEDTESTS_H

#include "Solvers/elasticsecondordersolver.h"
#include <fstream>
using namespace std;

class ElasticReducedTests
{
public:
    ElasticReducedTests();

    static void solveBartonTest1(int cellCount, int reinitialisationFrequency);

    static void solveZhangTest(int cellCount, int reinitialisationFrequency);

    static void solve2DBartonTest1(int cellCount, int reinitialisationFrequency);

    static void solve2DZhangTest(int cellCount, int reinitialisationFrequency);

    static void outputSolution(vector<ElasticReducedStateVector> solution);
    static void outputSolution2D(vector<vector<ElasticReducedStateVector> > solution);
};

#endif // ELASTICREDUCEDTESTS_H
