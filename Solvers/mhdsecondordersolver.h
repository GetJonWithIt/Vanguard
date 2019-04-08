#ifndef MHDSECONDORDERSOLVER_H
#define MHDSECONDORDERSOLVER_H

#include "mhdfirstordersolver.h"
using namespace std;

class MHDSecondOrderSolver
{
public:
    MHDSecondOrderSolver();

    static vector<double> computeXSLICFlux(MHDStateVector leftLeftStateVector, MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDStateVector rightRightStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);

    static void computeSLICTimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                    int slopeLimiter, MHDMaterialParameters materialParameters);

    static vector<MHDStateVector> solve(vector<MHDStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                        int subcyclingIterations, MHDMaterialParameters materialParameters);
};

#endif // MHDSECONDORDERSOLVER_H
