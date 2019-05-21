#ifndef HPRSOLVERS_H
#define HPRSOLVERS_H

#include "secondordersolver.h"
#include "HPR/hpracoustictensor.h"
#include "HPR/Multiphysics/hprreducedacoustictensor.h"
using namespace std;

class HPRSolvers
{
public:
    HPRSolvers();

    static vector<HPRStateVector> insertBoundaryCells(vector<HPRStateVector> & currentCells, int boundarySize);
    static vector<HPRReducedStateVector> insertBoundaryCells(vector<HPRReducedStateVector> & currentCells, int boundarySize);

    static vector<vector<HPRStateVector> > insertBoundaryCells2D(vector<vector<HPRStateVector> > & currentCells, int boundarySize);
    static vector<vector<HPRReducedStateVector> > insertBoundaryCells2D(vector<vector<HPRReducedStateVector> > & currentCells, int boundarySize);

    static double computeMaximumWaveSpeed(vector<HPRStateVector> & currentCells, HPRMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed(vector<HPRReducedStateVector> & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static double computeMaximumWaveSpeed2D(vector<vector<HPRStateVector> > & currentCells, HPRMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed2D(vector<vector<HPRReducedStateVector> > & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static double computeStableTimeStep(vector<HPRStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        HPRMaterialParameters materialParameters);
    static double computeStableTimeStep(vector<HPRReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static double computeStableTimeStep2D(vector<vector<HPRStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                          HPRMaterialParameters materialParameters);
    static double computeStableTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static HPRStateVector evolveStateByHalfXTimeStep(HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, HPRMaterialParameters materialParameters);
    static HPRReducedStateVector evolveStateByHalfXTimeStep(HPRReducedStateVector leftStateVector, HPRReducedStateVector middleStateVector, HPRReducedStateVector rightStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, HPRMaterialParameters material1Parameters,
                                                            HPRMaterialParameters material2Parameters);

    static HPRStateVector evolveStateByHalfYTimeStep(HPRStateVector topStateVector, HPRStateVector middleStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, HPRMaterialParameters materialParameters);
    static HPRReducedStateVector evolveStateByHalfYTimeStep(HPRReducedStateVector topStateVector, HPRReducedStateVector middleStateVector, HPRReducedStateVector bottomStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, HPRMaterialParameters material1Parameters,
                                                            HPRMaterialParameters material2Parameters);

    static HPRStateVector evolveStateByFractionalXTimeStep(double stepFraction, HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector,
                                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static HPRReducedStateVector evolveStateByFractionalXTimeStep(double stepFraction, HPRReducedStateVector leftStateVector, HPRReducedStateVector middleStateVector,
                                                                  HPRReducedStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                  HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static HPRStateVector evolveStateByFractionalYTimeStep(double stepFraction, HPRStateVector topStateVector, HPRStateVector middleStateVector, HPRStateVector bottomStateVector, double cellSpacing,
                                                           double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static HPRReducedStateVector evolveStateByFractionalYTimeStep(double stepFraction, HPRReducedStateVector topStateVector, HPRReducedStateVector middleStateVector,
                                                                  HPRReducedStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                  HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static HPRStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                    HPRMaterialParameters materialParameters);
    static HPRReducedStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static HPRStateVector evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution, HPRMaterialParameters materialParameters);
    static HPRReducedStateVector evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                 HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
};

#endif // HPRSOLVERS_H
