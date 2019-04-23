#ifndef HPRSOLVERS_H
#define HPRSOLVERS_H

#include "secondordersolver.h"
#include "HPR/hpracoustictensor.h"
using namespace std;

class HPRSolvers
{
public:
    HPRSolvers();

    static vector<HPRStateVector> insertBoundaryCells(vector<HPRStateVector> & currentCells, int boundarySize);

    static vector<vector<HPRStateVector> > insertBoundaryCells2D(vector<vector<HPRStateVector> > & currentCells, int boundarySize);

    static double computeMaximumWaveSpeed(vector<HPRStateVector> & currentCells, HPRMaterialParameters materialParameters);

    static double computeMaximumWaveSpeed2D(vector<vector<HPRStateVector> > & currentCells, HPRMaterialParameters materialParameters);

    static double computeStableTimeStep(vector<HPRStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        HPRMaterialParameters materialParameters);

    static double computeStableTimeStep2D(vector<vector<HPRStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                          HPRMaterialParameters materialParameters);

    static HPRStateVector evolveStateByHalfXTimeStep(HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, HPRMaterialParameters materialParameters);

    static HPRStateVector evolveStateByHalfYTimeStep(HPRStateVector topStateVector, HPRStateVector middleStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, HPRMaterialParameters materialParameters);

    static HPRStateVector evolveStateByFractionalXTimeStep(double stepFraction, HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector,
                                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);

    static HPRStateVector evolveStateByFractionalYTimeStep(double stepFraction, HPRStateVector topStateVector, HPRStateVector middleStateVector, HPRStateVector bottomStateVector, double cellSpacing,
                                                           double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);

    static HPRStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                    HPRMaterialParameters materialParameters);

    static HPRStateVector evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution, HPRMaterialParameters materialParameters);
};

#endif // HPRSOLVERS_H
