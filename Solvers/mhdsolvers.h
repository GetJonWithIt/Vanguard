#ifndef MHDSOLVERS_H
#define MHDSOLVERS_H

#include "secondordersolver.h"
using namespace std;

class MHDSolvers
{
public:
    MHDSolvers();

    static vector<MHDStateVector> insertBoundaryCells(vector<MHDStateVector> & currentCells, int boundarySize);

    static double computeMaximumWaveSpeed(vector<MHDStateVector> & currentCells, MHDMaterialParameters materialParameters);

    static double computeStableTimeStep(vector<MHDStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        MHDMaterialParameters materialParameters);

    static MHDStateVector evolveStateByHalfXTimeStep(MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, MHDMaterialParameters materialParameters);

    static MHDStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                    MHDMaterialParameters materialParameters);
};

#endif // MHDSOLVERS_H
