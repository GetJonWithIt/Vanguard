#ifndef EULERTESTS_H
#define EULERTESTS_H

#include "Solvers/secondordersolver.h"
#include <fstream>
using namespace std;

class EulerTests
{
public:
    EulerTests();

    static void solveToroTest1(int cellCount);
    static void solveToroTest2(int cellCount);
    static void solveToroTest3(int cellCount);
    static void solveToroTest4(int cellCount);
    static void solveToroTest5(int cellCount);

    static void solve2DToroTest1(int cellCount);

    static void outputSolution(vector<EulerStateVector> solution);
    static void outputSolution2D(vector<vector<EulerStateVector> > solution);
};

#endif // EULERTESTS_H
