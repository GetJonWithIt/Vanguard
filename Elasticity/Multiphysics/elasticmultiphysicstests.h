#ifndef ELASTICMULTIPHYSICSTESTS_H
#define ELASTICMULTIPHYSICSTESTS_H

#include "Solvers/elasticsecondordersolver.h"
#include <fstream>
using namespace std;

class ElasticMultiphysicsTests
{
public:
    ElasticMultiphysicsTests();

    static void solveBartonTest1(int cellCount, int reinitialisationFrequency);
    static void solveBartonTest2(int cellCount, int reinitialisationFrequency);

    static void solveZhangTest(int cellCount, int reinitialisationFrequency);

    static void outputSolution(vector<ElasticMultiphysicsStateVector> solution);
};

#endif // ELASTICMULTIPHYSICSTESTS_H
