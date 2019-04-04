#ifndef ELASTICTESTS_H
#define ELASTICTESTS_H

#include "Solvers/secondordersolver.h"
#include <fstream>
using namespace std;

class ElasticTests
{
public:
    ElasticTests();

    static void solveBartonTest1(int cellCount);
    static void solveBartonTest2(int cellCount);

    static void solve2DBartonTest1(int cellCount);
    static void solve2DBartonTest2(int cellCount);

    static void outputSolution(vector<ElasticStateVector> solution);
    static void outputSolution2D(vector<vector<ElasticStateVector> > solution);
};

#endif // ELASTICTESTS_H
