#ifndef MHDSOLVERS_H
#define MHDSOLVERS_H

#include "secondordersolver.h"
using namespace std;

class MHDSolvers
{
public:
    MHDSolvers();

    static vector<MHDStateVector> insertBoundaryCells(vector<MHDStateVector> & currentCells, int boundarySize);
    static vector<MHDReducedStateVector> insertBoundaryCells(vector<MHDReducedStateVector> & currentCells, int boundarySize);

    static vector<vector<MHDStateVector> > insertBoundaryCells2D(vector<vector<MHDStateVector> > & currentCells, int boundarySize);
    static vector<vector<MHDReducedStateVector> > insertBoundaryCells2D(vector<vector<MHDReducedStateVector> > & currentCells, int boundarySize);

    static double computeMaximumWaveSpeed(vector<MHDStateVector> & currentCells, MHDMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed(vector<MHDReducedStateVector> & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeMaximumWaveSpeed2D(vector<vector<MHDStateVector> > & currentCells, MHDMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed2D(vector<vector<MHDReducedStateVector> > & currentCells, MHDMaterialParameters material1Parameters,
                                            MHDMaterialParameters material2Parameters);

    static double computeStableTimeStep(vector<MHDStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        MHDMaterialParameters materialParameters);
    static double computeStableTimeStep(vector<MHDReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeStableTimeStep2D(vector<vector<MHDStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, MHDMaterialParameters materialParameters);
    static double computeStableTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByHalfXTimeStep(MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, MHDMaterialParameters materialParameters);
    static MHDReducedStateVector evolveStateByHalfXTimeStep(MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector rightStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, MHDMaterialParameters material1Parameters,
                                                            MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByHalfYTimeStep(MHDStateVector topStateVector, MHDStateVector middleStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, MHDMaterialParameters materialParameters);
    static MHDReducedStateVector evolveStateByHalfYTimeStep(MHDReducedStateVector topStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector bottomStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, MHDMaterialParameters material1Parameters,
                                                            MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByFractionalXTimeStep(double stepFraction, MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector,
                                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static MHDReducedStateVector evolveStateByFractionalXTimeStep(double stepFraction, MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector,
                                                                  MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                  MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByFractionalYTimeStep(double stepFraction, MHDStateVector topStateVector, MHDStateVector middleStateVector, MHDStateVector bottomStateVector,
                                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static MHDReducedStateVector evolveStateByFractionalYTimeStep(double stepFraction, MHDReducedStateVector topStateVector, MHDReducedStateVector middleStateVector,
                                                                  MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                  MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                    MHDMaterialParameters materialParameters);
    static MHDReducedStateVector evolveStateByHalfTimeStepReduced(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                  MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution, MHDMaterialParameters materialParameters);
    static MHDReducedStateVector evolveStateByFractionalTimeStepReduced(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                        MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDSOLVERS_H
